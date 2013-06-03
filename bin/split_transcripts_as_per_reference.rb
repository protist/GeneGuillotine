#!/usr/bin/env ruby
#
# This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
# This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# This script reads in a gtf file, and splits it into non-overlapping genes,
#   with one gene per gene according to a second (reference) gff. This might be
#   useful if a downstream program (e.g. DEXSeq) requires each gene to be
#   separate and not overlapping with its neighbours.
#
# The position of the split is decided by priority:
#   1) If there is an existing split in any transcript isoform within the gtf,
#        then use that.
#   2) Otherwise, use the local minimum for depth of coverage (from a pileup
#        file created by samtools mpileup).
#   2b) If there are multiple minima, then use the one closest to the middle of
#         the gap.
#
# Output file:
#   gtf with gene IDs from reference gff (N.B. DEXSeq does not show position of
#   transcript, only gene ID).
#
# Notes on pileup:
#   To create a master pileup from several bam files, first merge bam files with
#   `samtools merge merged.bam 1.bam 2.bam etc.`, then create pileup with
#   `samtools mpileup merged.bam > merged.pileup`

require 'optparse'
options = {}
OptionParser.new do |opts|
  opts.banner='This script reads in a gtf file, and splits it into non-'\
              'overlapping genes, with one gene per gene according to a second '\
              '(reference) gff. This might be useful if a downstream program '\
              '(e.g. DEXSeq) requires each gene to be separate and not '\
              "overlapping with its neighbours.\nUsage:"
  opts.on_tail('-h', '--help', 'Show this message') do
    puts opts; exit
  end
  opts.on('-v', '--[no-]verbose', 'Run verbosely') do |v|
    options[:verbose] = v
  end
  opts.on('-i', '--input GTF_FILE',
          'Primary GTF_FILE (e.g. from Cuffmerge) is required.') do |i|
    options[:mygtf_path] = i
  end
  opts.on('-g', '--ref_gff GFF_FILE',
          'Reference GFF_FILE (e.g. from eupathdb) is required.') do |g|
    options[:refgff_path] = g
  end
  opts.on('-p', '--pileup PILEUP_FILE',
          'PILEUP_FILE (from samtools mpileup) is required.') do |p|
    options[:pileup_path] = p
  end
  opts.on('-o', '--output OUTPUT_FILE',
          'OUTPUT_FILE is required.') do |o|
    options[:output_path] = o
  end
end.parse!

# Mandatory "options"
raise OptionParser::MissingArgument if options[:mygtf_path].nil? ||
    options[:refgff_path].nil? ||
    options[:pileup_path].nil? ||
    options[:output_path].nil?

# Ruby 1.9.2 is required for ordered hash. Otherwise, I could hash.sort every time. But I won't.
min_release = '1.9.2'
ruby_release = RUBY_VERSION
if ruby_release < min_release
  abort("This script requires Ruby version #{min_release} or later. You are running #{ruby_release}.")
end

################################################################################
### Define string modifiers when transcripts aren't on one gene. These can
###   mostly be modified here, except were indicated.
class Symbol
# Define the gene ID to use if the transcript does not overlap any genes.
  # Transcript is intergenic.
  def between(downstream_gene)
    (self.to_s + '_to_' + downstream_gene.to_s).to_sym
  end
  # Transcript is after all annotated reference genes on this contig.
  def end
    ('after_last_' + self.to_s).to_sym
  end
  # Transcript is before all annotated reference genes on this contig.
  def begin
    ('before_first_' + self.to_s).to_sym
  end

# Define the transcript ID to use if the transcript needs to be split.
  def multiple
    second_suffix = ':2'
    suffix_regex = /:[0-9]*$/
    if !(self.to_s =~ suffix_regex)
      (self.to_s + second_suffix).to_sym
    else
      self.to_s.next.to_sym
    end
  end
end

################################################################################
### Read primary gtf file from user's experimental data
# Define UserTranscripts class
class UserTranscripts
  # A UserTranscripts object stores information from a cuffmerge output, specifically
  # the chromosome, transcript ID, start and stop coordinates of each exon, and
  # "other" information (e.g. strand, oId, tss_id)
  # This object will also store statistics on how many transcripts overlap with
  #   none, one or x genes, or are on a reference contigs with no genes.
  def initialize(transcripts_by_chromosome = {}, overlap_stats = {})
    @transcripts_by_chromosome = transcripts_by_chromosome
    @overlap_stats = overlap_stats # e.g. {NA=>1, 0=>2, 1=>20, 2=>5, 3=>2, 5=>1}
                                   # where NA is for no genes, and 0 encompasses
                                   # both intergenic and terminal.
    @transcripts_to_split = {}
    @previously_split_transcripts = {} # {:parent_transcript_id => :transcript_id:last#}
    @split_coord_from_phase_one = {}
  end

  # Create a new chromosome and gene id if necessary, and replace start and stop
  # coordinates if they increase the boundaries.
  # {:chromosome => {:transcript_id => {
  #   :coords => [[start, end], [start, end]…],
  #   :other => ["+/-/.", "\"oId ... to end"]}}}
  # Then later add :gene_id for each transcript.
  def write_transcript(chromosome, transcript_id, start, stop, strand, notes)
    # Possible speedup: rather than checking for existence of chromosome and
    #   transcript, etc., just store value into a persistent variable and check
    #   for change (i.e. create it the first time, then just write subsequent times.)
    # Create new chromosome hash if it doesn't exist.
    @transcripts_by_chromosome[chromosome.to_sym] ||= {}
    # Create new transcript hash if it doesn't exist.
    @transcripts_by_chromosome[chromosome.to_sym][transcript_id.to_sym] ||= {}
    # Create new coord array if it doesn't exist.
    @transcripts_by_chromosome[chromosome.to_sym][transcript_id.to_sym][:coords] ||= []
    # Append coordinates to this array.
    @transcripts_by_chromosome[chromosome.to_sym][transcript_id.to_sym][:coords].push([start, stop])
    # Create and write to new :other array if it doesn't exist.
    @transcripts_by_chromosome[chromosome.to_sym][transcript_id.to_sym][:other] ||= [strand, notes]
  end

  attr_reader(:transcripts_by_chromosome, :transcripts_to_split, :previously_split_transcripts) # for debugging

  # For each chromosome, sort genes by start coordinates.
  # TODO: sort coordinates?
  def sort!
    @transcripts_by_chromosome.each do |chromosome, transcripts_for_this_chromosome|
      sorted_chromosome = Hash[transcripts_for_this_chromosome.sort_by { |_, value| value[:coords].first.first }]
      if sorted_chromosome.keys == transcripts_for_this_chromosome.keys
        puts "#{Time.new}:   chromosome #{chromosome} was already ordered correctly."
      else
        puts "#{Time.new}:   chromosome #{chromosome} is now ordered correctly."
        @transcripts_by_chromosome[chromosome] = sorted_chromosome
      end
    end
  end

  # For a particular chromosome and transcript_id, return exon coords in array of arrays.
  def transcript_coords(chromosome, transcript_id)
    @transcripts_by_chromosome[chromosome][transcript_id][:coords]
  end

  # Return names of the chromosomes as an array.
  def chromosome_names
    @transcripts_by_chromosome.keys
  end

  # Return transcript_ids for a given chromosome, as an array.
  def transcript_ids(chromosome)
    @transcripts_by_chromosome[chromosome].keys
  end

  # Add to statistics about where transcripts lie.
  #   e.g. {:NA=>1, 0=>2, 1=>20, 2=>5, 3=>2, 5=>1},
  #   with :NA for no ref contig, and 0 for both intergenic and terminal.
  def add_event(num_genes_overlap)
    @overlap_stats[num_genes_overlap] ||= 0
    @overlap_stats[num_genes_overlap] += 1
  end

  # Output statistics about where the transcripts lie.
  #   e.g. {:NA=>1, 0=>2, 1=>20, 2=>5, 3=>2, 5=>1},
  #   with :NA for no ref contig, and 0 for both intergenic and terminal.
  def overlap_stats
    output_hash = @overlap_stats.select { |key, _| key==:NA}
    output_hash.merge(Hash[@overlap_stats.select {|key, _| key!=:NA}.
                               sort_by {|key, _| key }])
  end

  # Write gene id to a transcript.
  def write_gene_id(chromosome, transcript_id, gene_id)
    @transcripts_by_chromosome[chromosome][transcript_id][:gene_id] = gene_id
  end

  # Add to list of transcripts to split later. Each member of :upstream_gene_id
  #   corresponds with a member of :coord.
  def define_split(chromosome, transcript_id, split_coord, upstream_gene_id)
    @transcripts_to_split[chromosome] ||= {}
    @transcripts_to_split[chromosome][transcript_id] ||= {coords:[], upstream_gene_ids:[]}
    @transcripts_to_split[chromosome][transcript_id][:coords].push(split_coord)
    @transcripts_to_split[chromosome][transcript_id][:upstream_gene_ids].push(upstream_gene_id)
  end

  # Find next transcript_id to use when splitting transcripts, and record
  #   biggest used transcript_id.
  def advance_transcript_id(chromosome, base_transcript_id, just_used_transcript_id)
    @previously_split_transcripts[chromosome] ||= {}
    if @previously_split_transcripts[chromosome][base_transcript_id]
      previously_last_used_transcript_id = @previously_split_transcripts[chromosome][base_transcript_id]
      previously_first_unused_transcript_id = previously_last_used_transcript_id.multiple
    else
      previously_first_unused_transcript_id = just_used_transcript_id
    end
    if just_used_transcript_id == previously_first_unused_transcript_id # i.e. used the next available transcript id.
      @previously_split_transcripts[chromosome][base_transcript_id] = just_used_transcript_id
      just_used_transcript_id.multiple
    else # i.e. had used an existing transcript id.
      previously_first_unused_transcript_id
    end
  end

  # Record split (from phase one, to potentially reuse in phase two).
  def record_split_coord(chromosome, upstream_gene_id, coord)
    @split_coord_from_phase_one[chromosome] ||= {}
    @split_coord_from_phase_one[chromosome][upstream_gene_id] = coord
  end

  # Read split (from phase one), otherwise return nil.
  def read_phase_one_split_coords(chromosome, upstream_gene_id)
    if @split_coord_from_phase_one[chromosome]
      @split_coord_from_phase_one[chromosome][upstream_gene_id]
    else
      nil
    end
  end

  # Make the splits. Call this after the loop is completed, to prevent problems
  #   with changing order of hash, and having to unnecessarily read new
  #   positions from the end.
  # TODO: don't split if directly adjacent to CDS. Read in option to specify
  #   minimum distance to CDS, either as percent of the intron, or an absolute
  #   value.
  # If split_coord lies outside transcript, do nothing.
  # TODO: if you keep splitting a transcript until there is nothing left, then
  #   there will be an error. But this is surely implausible?
  def split!
    @transcripts_to_split.each do |chromosome, transcripts_to_split_by_chromosome|
      transcripts_to_split_by_chromosome.each do |parent_transcript_id, split_info|
        working_transcript = @transcripts_by_chromosome[chromosome].delete(parent_transcript_id)
        if !working_transcript[:base_transcript_id]
          working_transcript[:base_transcript_id] = parent_transcript_id
        end
        new_transcript_id = parent_transcript_id
        split_info[:coords].each_with_index do |split_coord, index|
          if split_coord.between?(working_transcript[:coords].first.first, \
              working_transcript[:coords].last.last) # split within transcript
            found_split = false
            new_transcript_coords = []
            while !found_split
              if split_coord.between?(working_transcript[:coords].first.first, \
                  working_transcript[:coords].first.last) # in exon
                found_split = true
                # Add this exon if its length is non-zero.
                if !(split_coord == working_transcript[:coords].first.first)
                  new_transcript_coords.push([working_transcript[:coords].first.first, (split_coord - 1)])
                end
                # How much of the working transcript do we retain?
                if split_coord == working_transcript[:coords].first.last
                  working_transcript[:coords].shift
                else
                  working_transcript[:coords][0][0] = split_coord + 1
                end
              else
                # Check next intron. (Error following the last exon, but won't occur.)
                if split_coord.between?((working_transcript[:coords].first.last + 1), \
                    (working_transcript[:coords][1].first - 1))
                  found_split = true
                end
                # Not in this exon, so move these coords to transcript.
                new_transcript_coords.push(working_transcript[:coords].shift)
              end
            end
            # Write this just-split transcript if non-zero.
            if !(new_transcript_coords == [])
              @transcripts_by_chromosome[chromosome][new_transcript_id] = \
                  {coords: new_transcript_coords, other: working_transcript[:other], base_transcript_id: working_transcript[:base_transcript_id]}
              self.write_gene_id(chromosome, new_transcript_id, split_info[:upstream_gene_ids][index])
              new_transcript_id = self.advance_transcript_id(chromosome, working_transcript[:base_transcript_id], new_transcript_id)
            end
          elsif split_coord > working_transcript[:coords].last.last # hack to correct gene ID
            working_transcript[:gene_id] = split_info[:upstream_gene_ids][index]
          end
        end
        # Write last transcript. N.B. gene ID already set from parent transcript.
        # TODO: if non-zero?
        @transcripts_by_chromosome[chromosome][new_transcript_id] = working_transcript
        self.advance_transcript_id(chromosome, working_transcript[:base_transcript_id], new_transcript_id)
      end
    end
    @transcripts_to_split = {}
  end

  # Return transcript ids that have a given gene id. If no transcripts for that
  #   gene (or chromosome), return nil. Also return the minimum and maximum
  #   coordinates for these transcripts. Return as [array of transcripts,
  #   [min/max coords]]
  def transcripts_and_coords_union_for_gene(chromosome, gene_id)
    if @transcripts_by_chromosome[chromosome]
      transcripts_for_this_gene = @transcripts_by_chromosome[chromosome].select do |_, transcript|
        transcript[:gene_id] == gene_id
      end
      if transcripts_for_this_gene == []
        nil
      else
        union_coords = [Float::INFINITY , -Float::INFINITY]
        transcripts_for_this_gene.each do |_, transcripts|
          if transcripts[:coords].first.first < union_coords.first
            union_coords[0] = transcripts[:coords].first.first
          end
          if union_coords.last < transcripts[:coords].last.last
            union_coords[1] = transcripts[:coords].last.last
          end
        end
        [transcripts_for_this_gene, union_coords]
      end
    else
      nil
    end
  end

  # Write to a file in the same format as the input gtf.
  #   I won't worry about writing exon number nor gene name for now.
  def write_to_file(output_path)
    File.open(output_path, 'w') do |output_file|
      @transcripts_by_chromosome.each do |chromosome_name, transcripts|
        transcripts.each do |transcript_id, transcript_info|
          other = transcript_info[:other]
          transcript_info[:coords].each do |exon_coords|
            output_line = [chromosome_name, 'Cufflinks', 'exon', \
                exon_coords.first, exon_coords.last, '.', other[0], '.', \
                "gene_id \"#{transcript_info[:gene_id]}\"; "\
                "transcript_id \"#{transcript_id}\"; #{other[1]}"]
            output_file.puts(output_line.join("\t"))
          end
        end
      end
    end
  end
end

# Read transcripts from cuffmerge out
#   chr <Cufflinks> <exon> start end <.> +/-/. <.> gene_id "XLOC_000384";
#     transcript_id "TCONS_00001640"; exon_number "3"; oId "CUFF.506.1"; tss_id "TSS678";
#     occasionally, ends with oId "CUFF.1863.2"; contained_in "TCONS_00005635"; tss_id "TSS2590";
#     Also occasionally has nearest_ref and class_code after oId.
#     if strand == "-", exons are still in increasing order.
#   everything in <> is consistent across the whole file, so don't bother storing it.
#   chr and transcript_id will be stored as a key.
#   gene_id will be replaced from id from reference gff.
#   exon_number will be implied by the order in the array.
#   +/-/., oId and tss_id are unnecessary for DEXSeq (but let's hang on to them
#     for now, unless this has a large effect on efficiency).
puts "#{Time.new}: parsing primary gtf file."
transcripts = UserTranscripts.new
File.open(options[:mygtf_path]).each do |line| # There are no header lines for cuffmerge out.
                                               # These are the parts of each line that we need.
  splitline = line.chomp.split("\t")
  if splitline !=[]
    transcripts.write_transcript(splitline[0], \
        /transcript_id "([^"]*)"/.match(splitline[8])[1], \
        splitline[3].to_i, splitline[4].to_i, splitline[6], \
        /oId ".*$/.match(splitline[8])[0])
  end
end

# Check that file is ordered
puts "#{Time.new}: checking order of primary gtf file."
transcripts.sort!

################################################################################
### Read reference GFF file (e.g. from eupathdb).
# Define ReferenceGFF class
class ReferenceGFF
  # A ReferenceGFF object stores reference gene models, specifically the
  # chromosome, gene ID, and terminal start and stop coordinates (e.g. of the CDS).
  def initialize(genes_by_chromosome = {})
    @genes_by_chromosome = genes_by_chromosome
  end

  attr_reader(:genes_by_chromosome) # for debugging

  # Create a new chromosome and gene id if necessary, and replace start and stop
  # coordinates if they increase the boundaries.
  # {:chromosome => {:geneID => [start, stop]} }
  def write_gene(chromosome, gene_id, start, stop)
    # Create new chromosome hash if it doesn't exist.
    @genes_by_chromosome[chromosome.to_sym] ||= {}
    # Create new gene_id hash if it doesn't exist.
    @genes_by_chromosome[chromosome.to_sym][gene_id.to_sym] ||= [start, stop]
    # Compare coordinates to existing ones (redundancy if you've just created
    #   the gene, but I'm not sure if it's still quicker to use ||= )
    if start < @genes_by_chromosome[chromosome.to_sym][gene_id.to_sym].first
      @genes_by_chromosome[chromosome.to_sym][gene_id.to_sym][0] = start
    end
    if stop > @genes_by_chromosome[chromosome.to_sym][gene_id.to_sym].last
      @genes_by_chromosome[chromosome.to_sym][gene_id.to_sym][1] = stop
    end
  end

  # For each chromosome, sort genes by start coordinates.
  def sort!
    @genes_by_chromosome.each do |chromosome, genes_for_this_chromosome|
      sorted_chromosome = Hash[genes_for_this_chromosome.sort_by { |_, coords| coords[0] }]
      if sorted_chromosome.keys == genes_for_this_chromosome.keys
        puts "#{Time.new}:   chromosome #{chromosome} was sorted correctly."
      else
        puts "#{Time.new}:   WARNING! Chromosome #{chromosome} was not sorted correctly (but now it is)."
        @genes_by_chromosome[chromosome] = sorted_chromosome
      end
    end
  end

  # Check to see if adjacent genes overlap (or touch).
  #   Originally, I attempted to fix this by splitting adjacent genes at the
  #   midpoint. There are a few issues with this. Firstly, this code will check
  #   for the presence of overlap, but may miss some multiple cases (e.g. if
  #   three genes overlap). Also, sometimes a quick fix of splitting at the
  #   midpoint is inappropriate. It's better to notify the user that there are
  #   issues and let them manually fix it (e.g. by deleting one of the genes).
  # N.B. an "overlap" of
  # For toxodb files, there are only 3 and 14 overlaps for GT1 and ME49
  #   respectively anyway.
  def check_overlaps
    overlap_count = 0
    total_count = 0
    @genes_by_chromosome.each do |chromosome, genes_for_this_chromosome|
      prev_gene_id = nil
      genes_for_this_chromosome.each do |gene_id, coords|
        if prev_gene_id
        total_count += 1
          if genes_for_this_chromosome[prev_gene_id].last >= coords.first - 1
            overlap_count += 1
            puts "#{Time.new}:   ERROR! On #{chromosome}, genes "\
              "#{prev_gene_id.to_s} and #{gene_id.to_s} overlap by (at least) "\
              "#{genes_for_this_chromosome[prev_gene_id].last - coords.first + 1}"\
              ' bp.'
          end
        end
        prev_gene_id = gene_id
      end
    end
    if overlap_count > 0
      abort("#{Time.new}:   ERROR! #{overlap_count} overlaps reported.")
    end
  end

  # Return names of the chromosomes as an array.
  def chromosome_names
    @genes_by_chromosome.keys
  end

  # Return transcript_ids for a given chromosome, as an array.
  def gene_ids(chromosome)
    @genes_by_chromosome[chromosome].keys
  end

  # Return the length of the chromosome.
  def chromosome_length(chromosome)
    if @genes_by_chromosome[chromosome]
      @genes_by_chromosome[chromosome].length
    else
      0
    end
  end

  # Return an array of the start and end coordinates of the gene, given
  # chromosome and gene_id
  def gene_coords(chromosome, gene_id)
    @genes_by_chromosome[chromosome][gene_id]
  end

  # Return an array of the start and end coordinates of the gene, given
  # chromosome and position in the ordered hash.
  def gene_coords_by_position(chromosome, position)
    # TODO: is it more efficient to store values and keys semi-permanently in separate arrays?
    @genes_by_chromosome[chromosome].values[position]
  end

  # Return the gene id, given chromosome and position in the ordered hash.
  def gene_id(chromosome, position)
    @genes_by_chromosome[chromosome].keys[position]
  end
end

# Since the UTRs from eupathdb are not all defined, best to be consistent and
#   define genes as being from first to last CDS (N.B. won't include rRNA).
# In field 9 -> Parent=rna_TGME49_203135-1 (N.B. all CDS are /-1$/).
# N.B. if strand == "-", arranged in reverse numerical order, but let's not make
#   either assumption here, just in case this file does not come from eupathdb.
# I think these files don't contain any alternatively-spliced genes. At least,
#   (for TGGT1 and TGME49) no entries have {$3 == "mRNA"} and contain "-2".
#   Even if they did, this script only considers the terminal exons of all CDS
#   sharing the same parent (without -1).
puts "#{Time.new}: parsing reference gff file."
refgff = ReferenceGFF.new
File.open(options[:refgff_path]).each do |line|
  splitline = line.split("\t")
  if splitline[2] == 'CDS' # Hence, ignore the header lines.
                           # These are the parts of the line that we need.
    refgff.write_gene(splitline[0], /Parent=rna_(.*)-1$/.match(splitline[8])[1], splitline[3].to_i, splitline[4].to_i)
  end
end
# TODO: also import rRNA from reference GFF

# Check that file is ordered.
puts "#{Time.new}: checking order of reference gff file."
refgff.sort!
puts "#{Time.new}: checking reference gff file for overlapping genes."
refgff.check_overlaps

################################################################################
### Read pileup file (from samtools mpileup)
# Define pileup class
class Pileup
  # A Pileup object stores relevant information from a pileup file, specifically
  #   chromosome, coordinates, and number of reads per coordinate.

  # A pileup file (from samtools mpileup)is tab-delimited, containing
  #   chromosome name, coordinate, reference base, the number of reads covering
  #   the site, read bases, and base qualities
  #   e.g. TGGT1_chrII  1  N  8 ^$C^$C^$C^$C^$C^$C^$C^$C  @@CC+@C@
  #   This is ordered by chromosome, then by coordinates. N.B. There are no
  #   entries for zero depth, but transcripts shouldn't contain these anyway.
  # The @pileup_by_chromosome has format as follows
  #   {:chromosome => {coordinate => depth, coordinate => depth, ...}}
  def initialize(pileup_path)
    @pileup_by_chromosome = {}
    File.open(pileup_path).each do |line|
      splitline = line.split "\t"
      input_chromosome = splitline[0].to_sym
      input_coordinate = splitline[1].to_i
      input_depth = splitline[3].to_i
      # Create new chromosome hash if it doesn't exist.
      @pileup_by_chromosome[input_chromosome] ||= {}
      # Add the coordinate and depth from this line of the pileup file.
      @pileup_by_chromosome[input_chromosome][input_coordinate] = input_depth
    end
  end

  # Find most-central minimum, given chromosome and coordinate range.
  # TODO: don't specify position directly adjacent to CDS. Add option to specify
  #   minimum distance to CDS, either as percent of the intron, or an absolute
  #   value.
  def minimum(chromosome, start, stop)
    max_coord_delta = ((stop - start)/2).to_i
    min_coord = nil
    min_depth = Float::INFINITY
    # Test from both extremities, moving towards the centre.
    (0..max_coord_delta).each do |coord_delta|
      [start + coord_delta, stop - coord_delta].each do |testing_coord|
        testing_depth = @pileup_by_chromosome[chromosome][testing_coord]
        if testing_depth <= min_depth
          min_coord = testing_coord
          min_depth = testing_depth
        end
      end
    end
    min_coord
  end
end

# Create pileup object
pileup = Pileup.new(options[:pileup_path])

################################################################################
### Find transcripts that overlap two genes; assign gene IDs to all.
# Split transcripts when they overlap two genes.
# There doesn't seem to be a way to change a key while retaining order.
#   Hence, the fixed array will no longer be sorted
# transcripts_for_this_chr["t3"] = transcripts_for_this_chr.delete("t1")

# Detect a problem if there a single gene straddles two in the reference.
#   Split the genes, and assign the correct gene ID to each.
# Walk through transcripts in order. Check from i=0 in refgff, until out of transcript.
#   Use the same i for the next gene. If the transcript does not fall in the
#   first j genes, increase i to (i + j) for the next transcript.
puts "#{Time.new}: fixing transcripts that overlap multiple reference genes."
transcripts.chromosome_names.each do |chromosome|
  puts "#{Time.new}:   analysing multiple genes per transcript for #{chromosome}."
  base_position_in_refgff_chr = 0
  length_of_refgff_chr = refgff.chromosome_length(chromosome)
  transcripts.transcript_ids(chromosome).each do |transcript_id|
    # For each transcript, find out into which (ref) genes the tss and tes fall.
    position_of_first_overlapping_gene = nil
    position_of_last_overlapping_gene = nil
    tss = transcripts.transcript_coords(chromosome, transcript_id).first.first
    tes = transcripts.transcript_coords(chromosome, transcript_id).last.last

    # Find first overlapping gene; technically, the first gene downstream of the
    #   tss, even if the transcript does not overlap any genes.
    # position_of_first_overlapping_gene = nil if transcript is downstream of all
    #   genes, or if there are no references genes for this chromosome.
    testing_position_in_refgff_chr = base_position_in_refgff_chr
    while !position_of_first_overlapping_gene && (testing_position_in_refgff_chr <= (length_of_refgff_chr - 1))
      testing_coords = refgff.gene_coords_by_position(chromosome, testing_position_in_refgff_chr)
      if tss > testing_coords.last # i.e. tss is downstream of gene, so reiterate.
        testing_position_in_refgff_chr += 1
        base_position_in_refgff_chr += 1
      else # i.e. tss is within gene, or upstream of gene.
        position_of_first_overlapping_gene = testing_position_in_refgff_chr
      end
    end

    # Find the last overlapping gene; technically, the last gene upstream of the
    #   tes, even if the transcript does not overlap any genes.
    # position_of_last_overlapping_gene = nil if transcript is downstream of all
    #   genes, or if there are no references genes for this chromosome.
    if testing_position_in_refgff_chr > 0
      testing_position_in_refgff_chr -=1 # for transcripts that overlap no genes.
    end
    while !position_of_last_overlapping_gene && (testing_position_in_refgff_chr <= (length_of_refgff_chr - 1))
      testing_coords = refgff.gene_coords_by_position(chromosome, testing_position_in_refgff_chr)
      if tes > testing_coords.last # i.e. tes is downstream of gene, but is it the last gene that overlaps?
        testing_position_in_refgff_chr += 1
      elsif tes >= testing_coords.first # i.e. tes is within gene
        position_of_last_overlapping_gene = testing_position_in_refgff_chr
      else # i.e. tes is upstream of gene
        position_of_last_overlapping_gene = testing_position_in_refgff_chr - 1
      end
    end
    if (!position_of_last_overlapping_gene) && (length_of_refgff_chr > 0) # i.e. at end of chromosome, but chromosome is not empty
      position_of_last_overlapping_gene = (length_of_refgff_chr - 1)
    end

    if options[:verbose]
      print "  #{transcript_id}: "
    end
    # What is the outcome of our tests?
    if !position_of_first_overlapping_gene
      if !position_of_last_overlapping_gene # no genes on this chromosome
        if options[:verbose]
          puts 'no genes on this reference contig'
        end
        transcripts.add_event(:NA)
        transcripts.write_gene_id(chromosome, transcript_id,
                                  :'No_genes_on_ref_contig')
      else # at the end of the chromosome
           # Could introduce this test earlier, and quickly mark all remaining
           #   transcripts identically, but there shouldn't be many, and it's
           #   not costly to iterate.
        if options[:verbose]
          puts "nearest_neighbour is the last gene, #{refgff.
              gene_id(chromosome, (length_of_refgff_chr - 1))}"
        end
        transcripts.add_event(0)
        transcripts.write_gene_id(chromosome, transcript_id, refgff.
            gene_id(chromosome, (length_of_refgff_chr - 1)).end)
      end
    elsif position_of_last_overlapping_gene == -1 # beginning of the chromosome
      if options[:verbose]
        puts "nearest_neighbour is the first gene, #{refgff.
            gene_id(chromosome, (0))}"
      end
      transcripts.add_event(0)
      transcripts.write_gene_id(chromosome, transcript_id, refgff.
          gene_id(chromosome, (0)).begin)
    elsif position_of_first_overlapping_gene == position_of_last_overlapping_gene
      if options[:verbose]
        puts "covers one gene, #{refgff.
            gene_id(chromosome, position_of_first_overlapping_gene)}"
      end
      transcripts.add_event(1)
      transcripts.write_gene_id(chromosome, transcript_id, refgff.
          gene_id(chromosome, position_of_first_overlapping_gene))
    elsif position_of_last_overlapping_gene < position_of_first_overlapping_gene # the transcript is wholly intergenic
      if options[:verbose]
        puts 'intergenic, between '\
            "#{refgff.gene_id(chromosome, position_of_last_overlapping_gene)} and "\
            "#{refgff.gene_id(chromosome, position_of_first_overlapping_gene)}"
      end
      transcripts.add_event(0)
      transcripts.write_gene_id(chromosome, transcript_id, \
          [refgff.gene_id(chromosome, position_of_last_overlapping_gene), \
          refgff.gene_id(chromosome, position_of_first_overlapping_gene)])
    else # Covers multiple genes. Find and record split positions.
      if options[:verbose]
        puts "covers #{position_of_last_overlapping_gene - \
          position_of_first_overlapping_gene + 1} genes from #{refgff.
            gene_id(chromosome, position_of_first_overlapping_gene)} to "\
          "#{refgff.gene_id(chromosome, position_of_last_overlapping_gene)}"
      end
      transcripts.add_event(position_of_last_overlapping_gene - position_of_first_overlapping_gene + 1)
      transcripts.write_gene_id(chromosome, transcript_id, refgff.
          gene_id(chromosome, position_of_last_overlapping_gene)) # Temporarily store this gene id for the entire transcript (easier when naming fragments later).
      (position_of_first_overlapping_gene..(position_of_last_overlapping_gene - 1)).
          each do |position_of_upstream_gene|
        split_coord = pileup.minimum(chromosome, (refgff.gene_coords_by_position(chromosome, position_of_upstream_gene).
            last + 1), (refgff.gene_coords_by_position(chromosome, position_of_upstream_gene + 1).first - 1))
        upstream_gene_id = refgff.gene_id(chromosome, position_of_upstream_gene)
        transcripts.define_split(chromosome, transcript_id, split_coord, upstream_gene_id)
        transcripts.record_split_coord(chromosome, upstream_gene_id, split_coord)
      end
    end
  end
end

transcripts.split!
puts "#{Time.new}: statistics for transcripts that overlap multiple reference genes."
transcripts.overlap_stats.each do |genes, count|
  if genes == :NA
    if count == 1
      puts '  1 transcript has no matching reference contig'
    else
      puts "  #{count} transcripts have no matching reference contig"
    end
  else
    if count == 1
      puts "  1 transcript overlaps with #{genes} genes"
    else
      puts "  #{count} transcripts overlap with #{genes} genes each"
    end
  end
end

################################################################################
### Find overlapping transcripts with different gene IDs
# TODO: what about ALAD-SPP? Only terminal exons? or not if they share the same tss? Getting complicated!
# Use same cut points as used previously. If the overlap still contains this
#   coordinate, they still fall in between the overlap, it makes sense to reuse
#   this. Hence, assume that the transcript overlap is less informative than the
#   positions of the reference genes. (I won't bother saving the previous split
#   points, unless recalculation is slow.) If the overlap does not contain the
#   split coordinate, find the minimal pileup position within the overlap.

puts "#{Time.new}: re-sorting transcripts."
transcripts.sort!

puts "#{Time.new}: fixing overlapping transcripts on adjacent genes."
refgff.chromosome_names.each do |chromosome|
  puts "#{Time.new}:   analysing adjacent overlapping transcripts for #{chromosome}."
  prev_gene_transcripts_and_coords = nil
  prev_gene_id = nil
  refgff.gene_ids(chromosome).each do |current_gene_id|
    current_gene_transcripts_and_coords = transcripts.
        transcripts_and_coords_union_for_gene(chromosome, current_gene_id)
    intergenic_transcripts_and_coords = transcripts.
        transcripts_and_coords_union_for_gene(chromosome, \
        [prev_gene_id, current_gene_id])
    # TODO: what if intergenic transcripts do not overlap with each other? (They should have been assigned unique gene ids.)
    split_coord = nil
    if prev_gene_id # Not the first iteration.
      if prev_gene_transcripts_and_coords && \
          current_gene_transcripts_and_coords
        # Is there overlap between the in-gene transcripts?
        if prev_gene_transcripts_and_coords.last.last >= \
          current_gene_transcripts_and_coords.last.first
          if options[:verbose]
            puts 'there is overlap between transcript(s) on genes ' \
                "#{prev_gene_id} and #{current_gene_id}"
          end
          # Check for a previous split in phase one between these genes.
          split_coord = transcripts.read_phase_one_split_coords(chromosome, prev_gene_id)
          # If this was not split on phase one, then find new coordinates
          if !split_coord
            # TODO: split_coord should be central in the intergenic region, not the overlap.
            # N.B. this finds a split coord within the overlap of transcripts
            #   between neighbouring genes, ignoring potential intergenic
            #   transcripts, which I consider "noisier".
            split_coord = pileup.minimum(chromosome, \
                current_gene_transcripts_and_coords.last.first, \
                prev_gene_transcripts_and_coords.last.last)
          end
        elsif intergenic_transcripts_and_coords && \
            prev_gene_transcripts_and_coords.last.last >= \
            intergenic_transcripts_and_coords.last.first && \
            intergenic_transcripts_and_coords.last.last >= \
            current_gene_transcripts_and_coords.last.first # intergenic overlaps both
          if options[:verbose]
            puts 'the intergenic transcript(s) overlap with transcript(s) on ' \
                "genes #{prev_gene_id} and #{current_gene_id}"
          end
          # Check for a previous split in phase one between these genes.
          split_coord = transcripts.read_phase_one_split_coords(chromosome, prev_gene_id)
          # If this was not split on phase one, then find new coordinates
          # The overlap union's coordinates are identical to the intergenic
          # transcripts'. I'm cutting as central as possible within this. There
          # may be an argument for not prioritising this, nor the whole inter-
          # genic region (as above), but it's unclear to me what would be the
          # best place to cut. Perhaps even outside the overlap in this case?
          # TODO: think about the best place to cut, as per this comment block.
          if !split_coord
            split_coord = pileup.minimum(chromosome, \
                intergenic_transcripts_and_coords.last.first, \
                intergenic_transcripts_and_coords.last.last)
          end
        end
        if split_coord
          # Cut prev_gene, intergenic, and current_gene transcripts.
          (prev_gene_transcripts_and_coords.first.keys + \
              intergenic_transcripts_and_coords.first.keys + \
              current_gene_transcripts_and_coords.first.keys).each do |transcript_id|
            transcripts.define_split(chromosome, transcript_id, split_coord, prev_gene_id)
            transcripts.write_gene_id(chromosome, transcript_id, current_gene_id)
          end
        end
      end
    end
    prev_gene_transcripts_and_coords = current_gene_transcripts_and_coords
    prev_gene_id = current_gene_id
  end
end

# TODO: if it's the terminal exon, just truncate, otherwise, split? N.B. different to behaviour of processing in first stage. Maybe don't throw out anything.
# TODO: perhaps just assign transcript ids at the end? In the meantime, just store in an array for the transcript id? Or just look for unique each time?

transcripts.split!

################################################################################
### Make gene name unique if transcripts do not overlap
# DEXSeq trusts geneIDs. Hence, it combines two genes if they have the same
#   geneID, regardless of where they are located.

# TODO: if multiple, overlapping transcripts on a single gene -> they should share a gene ID. We are trusting the refgff gene limits. Or should we totally trust the ref?
#   OTOH, if there are multiple non-overlapping transcripts on a single gene, we should break them up into -a, -b, etc.
#   similarly for intergenic transcripts.

#
# TODO: intergenic transcripts need to have their gene name numbered: e.g. :a, :b
# TODO: currently, intergenic transcripts have gene id as [:prev, :next].
#   Need to change this to string.between(second_string)
#   If intergenic transcripts overlap in-gene transcripts, just change their id to the relevant gene id?

puts 'SUMMARY
=======
'
puts transcripts.overlap_stats

transcripts.write_to_file(options[:output_path])