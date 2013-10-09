#!/usr/bin/env ruby
#
# Copyright 2013 Lee M. Yeoh (email: "plast-dot-id-dot-au" follows "github")
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
# The position of the split is determined by the genes in the reference gff. The
#   default is to constrain each transcript to the limits of the CDS.
#
# Transcripts that lie wholly within intergenic regions will be kept.
#
# Output file:
#   gtf with gene IDs from reference gff (N.B. DEXSeq does not show position of
#   transcript, only gene ID).

require 'optparse'
$options = {}
OptionParser.new do |opts|
  opts.banner='This script reads in a gtf file, and splits it into non-'\
              'overlapping genes, with one gene per gene according to a second '\
              '(reference) gff. This might be useful if a downstream program '\
              '(e.g. DEXSeq) requires each gene to be separate and not '\
              "overlapping with its neighbours.\nUsage:"
  opts.on_tail('-h', '--help', 'Show this message') do
    puts opts; exit
  end
  # Verbosity off (or == 0) is base-level information. Verbosity == 1 is
  #   chromosome-level. Verbosity == 2 is transcript/gene-level.
  opts.on('-v', '--[no-]verbose [OPT]', 'Run verbosely, optionally at level 2') do |v|
    $options[:verbosity] = (v || 1).to_i
  end
  opts.on('-i', '--input GTF_FILE',
        'Primary GTF_FILE (e.g. from Cuffmerge) is required.') do |i|
    $options[:mygtf_path] = i
  end
  opts.on('-g', '--ref_gff GFF_FILE',
        'Reference GFF_FILE (e.g. from eupathdb) is required.') do |g|
    $options[:refgff_path] = g
  end
  opts.on('-o', '--output OUTPUT_FILE',
        'OUTPUT_FILE is required.') do |o|
    $options[:output_path] = o
  end
  opts.on('-m', '--minimal_split', 'Only split transcripts that overlap ' \
      'multiple genes. Do not split when multiple adjacent transcripts ' \
      'overlap') do |m|
    $options[:minimal_split] = m
  end
end.parse!

# Mandatory "options"
raise OptionParser::MissingArgument if $options[:mygtf_path].nil? ||
    $options[:refgff_path].nil? ||
    $options[:output_path].nil?
$options[:verbosity] = 0 if !$options[:verbosity]

# TODO: allow option to allow an absolute extension on either side of the CDS
#   (e.g. 500 bp), while preventing adjacent defined UTRs from being too
#   close, perhaps by defining maximum UTR length? I'm not sure how useful this
#   would be in terms of power, it's slightly complicated to code (including
#   implementing all the options), and might result in more false positives.

# Require Ruby 1.9.2 for ordered hash. This is quicker than continual hash.sort.
min_release = '1.9.2'
ruby_release = RUBY_VERSION
if ruby_release < min_release
  abort("This script requires Ruby version #{min_release} or later. You are running #{ruby_release}.")
end

################################################################################
### Define string modifiers when transcripts aren't on one gene. These can
###   be modified.
class Symbol
# Define the gene ID to use if the transcript does not overlap any genes.
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

class Array
  # Transcript is intergenic. [:gene1, :gene2]
  def between
    if self.length == 2
      "[#{self.first}, #{self.last}]"
    else
      abort("#{self} is not an array of length two.")
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
  def initialize(transcripts_by_chromosome = {})
    @transcripts_by_chromosome = transcripts_by_chromosome
    @stats = {phase_one_overlaps:{}, phase_one_intergenic:{transcripts:0, \
        intergenics:0}, phase_two:{intergenics_adj_genes:0, \
        intergenics_adj_inter:0, transcripts:0}}
    @transcripts_to_split = {}
    @previously_split_transcripts = {} # {:parent_transcript_id => :transcript_id:last#}
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

  attr_reader(:transcripts_by_chromosome, :transcripts_to_split, :previously_split_transcripts, :stats)

  # For each chromosome, sort transcripts by start coordinates.
  # TODO: sort coordinates within transcripts? Not necessary, since the output
  #   from cuffmerge appears to always sort them already.
  def sort!
    @transcripts_by_chromosome.each do |chromosome, transcripts_for_this_chromosome|
      sorted_chromosome = Hash[transcripts_for_this_chromosome.sort_by { |_, value| value[:coords].first.first }]
      if sorted_chromosome.keys == transcripts_for_this_chromosome.keys
        if $options[:verbosity] >= 1
          puts "#{Time.new}:   chromosome #{chromosome} was already ordered correctly."
        end
      else
        if $options[:verbosity] >= 1
          puts "#{Time.new}:   chromosome #{chromosome} is now ordered correctly."
        end
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

  # Add to statistics about where transcripts lie and what has been processed.
  def add_event(category, data)
    case category
      when :phase_one_overlaps
        # data = the number of overlapping genes for a transcript
        #   e.g. {:NA=>1, 0=>2, 1=>20, 2=>5, 3=>2, 5=>1},
        #   with :NA for no ref contig, and 0 for both intergenic and terminal.
        @stats[:phase_one_overlaps][data] ||= 0
        @stats[:phase_one_overlaps][data] += 1
      when :phase_one_intergenic
        # data = [:intergenics, #] || [:transcripts, #]
        #   i.e. x number of intergenic transcripts removed over y intergenic
        #   regions. N.B. y refers to the amount of regions scrubbed, some might
        #   have been devoid of transcripts.
        @stats[:phase_one_intergenic][data.first] += data.last
      when :phase_two
        # data = [:intergenics_adj_genes, #] || [:intergenics_adj_inter] ||
        #   [:transcripts, #]
        #   i.e. in how many intergenic regions was overlap detected, and which
        #   kind was it? How many intergenic transcripts were deleted as a
        #   result?
        @stats[:phase_two][data.first] += data.last
      else
        abort("ERROR: add_event called with unknown category #{category}")
    end
  end

  # Write gene id to a transcript.
  def write_gene_id(chromosome, transcript_id, gene_id)
    @transcripts_by_chromosome[chromosome][transcript_id][:gene_id] = gene_id
  end

  # Add to list of transcripts to split later.
  #   @transcripts_to_split[chromosome][transcript_id] is an array containing
  #   [{coords:[a,b], upstream_gene_ids:[id_a], downstream_gene_id:[id_b]}, {…}]
  #   downstream_gene_id is only used for the first part of phase one.
  def define_split(chromosome, transcript_id, split_coords, upstream_gene_id, \
      downstream_gene_id)
    @transcripts_to_split[chromosome] ||= {}
    @transcripts_to_split[chromosome][transcript_id] ||= []
    @transcripts_to_split[chromosome][transcript_id].push\
        ({coords: split_coords, upstream_gene_id: upstream_gene_id, \
        downstream_gene_id: downstream_gene_id})
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

  # Make the splits. Call this after the loop is completed, to prevent problems
  #   with changing order of hash, and having to unnecessarily read new
  #   positions from the end.
  def split!
    @transcripts_to_split.each do |chromosome, transcripts_to_split_by_chromosome|
      transcripts_to_split_by_chromosome.each do |parent_transcript_id, split_events|
        working_transcript = @transcripts_by_chromosome[chromosome].delete(parent_transcript_id)
        if !working_transcript[:base_transcript_id]
          working_transcript[:base_transcript_id] = parent_transcript_id
        end
        new_transcript_id = parent_transcript_id
        more_to_analyse = true
        split_events.each do |split_event|
          if more_to_analyse
            # If transcript is totally upstream of split coords, don't split,
            #   fix gene_id, and stop more "splitting" and writing of gene_id.
            if working_transcript[:coords].last.last < split_event[:coords].first
              working_transcript[:gene_id] = split_event[:upstream_gene_id]
              more_to_analyse = false
            else

              # Firstly, see if upstream split position is in transcript.
              #   Then, split and create new transcript.
              if split_event[:coords].first.between? \
                  (working_transcript[:coords].first.first + 1), \
                  working_transcript[:coords].last.last
                new_transcript_coords = []
                found_split = false
                while !found_split
                  if split_event[:coords].first < working_transcript[:coords].first.first # upstream of exon
                    found_split = true
                  elsif split_event[:coords].first.between? \
                      working_transcript[:coords].first.first, \
                      working_transcript[:coords].first.last # in exon
                    found_split = true
                    # Add this exon if its length is non-zero.
                    if !(split_event[:coords].first == working_transcript[:coords].first.first)
                      new_transcript_coords.push([working_transcript[:coords].first.
                          first, (split_event[:coords].first - 1)])
                    end
                  else # in next intron or further downstream
                    # Check next intron. (Error following the last exon, but won't occur.)
                    if split_event[:coords].first.between? \
                        (working_transcript[:coords].first.last + 1), \
                        (working_transcript[:coords][1].first - 1)
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
                  self.write_gene_id(chromosome, new_transcript_id, split_event[:upstream_gene_id])
                  new_transcript_id = self.advance_transcript_id(chromosome, working_transcript[:base_transcript_id], new_transcript_id)
                end
              end

              # Secondly, see if downstream split position is in transcript.
              #   Then, trim this (possibly just-created) transcript.
              if working_transcript[:coords].last.last <= split_event[:coords].last
                more_to_analyse = false
                working_transcript[:coords] = []
              elsif split_event[:coords].last.between?(working_transcript[:coords].first.first, \
                  (working_transcript[:coords].last.last - 1))
                found_split = false
                while !found_split
                  if split_event[:coords].last < working_transcript[:coords].first.first # upstream of exon
                    found_split = true
                  elsif split_event[:coords].last.between?(working_transcript[:coords].first.first, \
                      working_transcript[:coords].first.last) # in exon
                    found_split = true
                    # Either remove exon or trim it if necessary.
                    if split_event[:coords].last == working_transcript[:coords].first.last
                      working_transcript[:coords].shift
                    else
                      working_transcript[:coords][0][0] = split_event[:coords].last + 1
                    end
                  else # In next intron or further downstream. Remove exon.
                    working_transcript[:coords].shift
                  end
                end
              end
            end
          end
        end
        # Write last working transcript if non-zero.
        #   gene_id will have been set from parent transcript.
        if !(working_transcript[:coords] == [])
          @transcripts_by_chromosome[chromosome][new_transcript_id] = working_transcript
          self.advance_transcript_id(chromosome, working_transcript[:base_transcript_id], new_transcript_id)
        end
      end
    end
    @transcripts_to_split = {}
  end

  # Return transcripts that have a given gene id as an array. If no transcripts
  #   for that gene (or chromosome), return empty array.
  def transcript_ids_for_gene(chromosome, gene_id)
    if @transcripts_by_chromosome[chromosome]
      @transcripts_by_chromosome[chromosome].select do |_, transcript_info|
        transcript_info[:gene_id] == gene_id
      end.collect do |transcript_id, _|
        transcript_id
      end
    else
      []
    end
  end

  # Return transcripts that have a given gene id. If no transcripts for that
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

  # Delete transcripts with a particular gene_id. Only used for intergenic
  #   transcripts. Returns the number of transcripts deleted.
  def delete_transcripts_for_gene(chromosome, gene_id)
    if @transcripts_by_chromosome[chromosome]
      initial_length = @transcripts_by_chromosome[chromosome].length
      @transcripts_by_chromosome[chromosome].delete_if do |_, transcript_info|
        transcript_info[:gene_id] == gene_id
      end
      initial_length - @transcripts_by_chromosome[chromosome].length
    else
      0
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
            gene_id = transcript_info[:gene_id]
            gene_id = gene_id.between if gene_id.class == Array
            output_line = [chromosome_name, 'Cufflinks', 'exon', \
                exon_coords.first, exon_coords.last, '.', other[0], '.', \
                "gene_id \"#{gene_id}\"; "\
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
File.open($options[:mygtf_path]).each do |line| # There are no header lines for cuffmerge out.
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
        if $options[:verbosity] >= 1
          puts "#{Time.new}:   chromosome #{chromosome} was sorted correctly."
        end
      else
        if $options[:verbosity] >= 1
          puts "#{Time.new}:   WARNING! Chromosome #{chromosome} was not sorted correctly (but now it is)."
        end
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
    @genes_by_chromosome[chromosome].values[position]
  end

  # Return the gene id, given chromosome and position in the ordered hash.
  def gene_id(chromosome, position)
    @genes_by_chromosome[chromosome].keys[position]
  end
end

# Since the UTRs from eupathdb are not all defined, best to be consistent and
#   define genes as being from first to last CDS (except for tRNA and rRNA).
# In field 9 -> Parent=rna_TGME49_203135-1 (N.B. all CDS are /-1$/).
# N.B. if strand == "-", arranged in reverse numerical order, but let's not make
#   either assumption here, just in case this file does not come from eupathdb.
# I think these files don't contain any alternatively-spliced genes. At least,
#   (for TGGT1 and TGME49) no entries have {$3 == "mRNA"} and contain "-2".
#   Even if they did, this script only considers the terminal exons of all CDS
#   sharing the same parent (without -1).
# GFF header lines won't have splitline[2].
puts "#{Time.new}: parsing reference gff file."
refgff = ReferenceGFF.new
File.open($options[:refgff_path]).each do |line|
  splitline = line.split("\t")
  skip = false
  if splitline[2] == 'CDS'
    matchdata = /Parent=rna_(TG[^_]{2,4}_\d*[A-Z]?)-1(;|$)/.match(splitline[8])
  elsif splitline[2] == 'tRNA' || splitline[2] == 'rRNA'
    matchdata = /Parent=(TG[^_]{2,4}_\d*)(;|$)/.match(splitline[8])
  else
    skip = true
  end
  if !skip
    if matchdata
      refgff.write_gene(splitline[0], matchdata[1], splitline[3].to_i, splitline[4].to_i)
    else
      abort('ERROR: the following line does not contain a gene_id in the ' \
          "expected format\n#{line}")
    end
  end
end

# Check that file is ordered.
puts "#{Time.new}: checking order of reference gff file."
refgff.sort!
puts "#{Time.new}: checking reference gff file for overlapping genes."
refgff.check_overlaps

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
  if $options[:verbosity] >= 1
    puts "#{Time.new}:   identifying multiple genes per transcript for #{chromosome}."
  end
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

    if $options[:verbosity] >= 2
      print "  #{transcript_id}: "
    end
    # What is the outcome of our tests?
    if !position_of_first_overlapping_gene
      if !position_of_last_overlapping_gene # no genes on this chromosome
        if $options[:verbosity] >= 2
          puts 'no genes on this reference contig'
        end
        transcripts.add_event(:phase_one_overlaps, :NA)
        transcripts.write_gene_id(chromosome, transcript_id,
                                  :'No_genes_on_ref_contig')
      else # at the end of the chromosome
           # Could introduce this test earlier, and quickly mark all remaining
           #   transcripts identically, but there shouldn't be many, and it's
           #   not costly to iterate.
        if $options[:verbosity] >= 2
          puts "nearest_neighbour is the last gene, #{refgff.
              gene_id(chromosome, (length_of_refgff_chr - 1))}"
        end
        transcripts.add_event(:phase_one_overlaps, 0)
        transcripts.write_gene_id(chromosome, transcript_id, refgff.
            gene_id(chromosome, (length_of_refgff_chr - 1)).end)
      end
    elsif position_of_last_overlapping_gene == -1 # beginning of the chromosome
      if $options[:verbosity] >= 2
        puts "nearest_neighbour is the first gene, #{refgff.
            gene_id(chromosome, (0))}"
      end
      transcripts.add_event(:phase_one_overlaps, 0)
      transcripts.write_gene_id(chromosome, transcript_id, refgff.
          gene_id(chromosome, (0)).begin)
    elsif position_of_first_overlapping_gene == position_of_last_overlapping_gene
      if $options[:verbosity] >= 2
        puts "covers one gene, #{refgff.
            gene_id(chromosome, position_of_first_overlapping_gene)}"
      end
      transcripts.add_event(:phase_one_overlaps, 1)
      transcripts.write_gene_id(chromosome, transcript_id, refgff.
          gene_id(chromosome, position_of_first_overlapping_gene))
    elsif position_of_last_overlapping_gene < position_of_first_overlapping_gene # the transcript is wholly intergenic
      if $options[:verbosity] >= 2
        puts 'intergenic, between '\
            "#{refgff.gene_id(chromosome, position_of_last_overlapping_gene)} and "\
            "#{refgff.gene_id(chromosome, position_of_first_overlapping_gene)}"
      end
      transcripts.add_event(:phase_one_overlaps, 0)
      transcripts.write_gene_id(chromosome, transcript_id, \
          [refgff.gene_id(chromosome, position_of_last_overlapping_gene), \
          refgff.gene_id(chromosome, position_of_first_overlapping_gene)])
    else # Covers multiple genes. Find and record split positions.
      if $options[:verbosity] >= 2
        puts "covers #{position_of_last_overlapping_gene - \
          position_of_first_overlapping_gene + 1} genes from #{refgff.
            gene_id(chromosome, position_of_first_overlapping_gene)} to "\
          "#{refgff.gene_id(chromosome, position_of_last_overlapping_gene)}"
      end
      transcripts.add_event(:phase_one_overlaps, \
          position_of_last_overlapping_gene - \
          position_of_first_overlapping_gene + 1)
      transcripts.write_gene_id(chromosome, transcript_id, refgff.
          gene_id(chromosome, position_of_last_overlapping_gene)) # Temporarily store this gene id for the entire transcript (easier when naming fragments later).
      (position_of_first_overlapping_gene..(position_of_last_overlapping_gene - 1)).
          each do |position_of_upstream_gene|
        split_coords = [(refgff.gene_coords_by_position(chromosome,
            position_of_upstream_gene).last + 1), (refgff.gene_coords_by_position(
            chromosome, position_of_upstream_gene + 1).first - 1)]
        upstream_gene_id = refgff.gene_id(chromosome, position_of_upstream_gene)
        downstream_gene_id = refgff.gene_id(chromosome, position_of_upstream_gene + 1)
        transcripts.define_split(chromosome, transcript_id, split_coords, \
            upstream_gene_id, downstream_gene_id)
      end
    end
  end
end

# For intergenic regions that will be excised because of the presence of
#   transcripts covering multiple genes, also remove intergenic transcripts and
#   extensions from other transcripts into this region.
# Don't integrate this identification into split!, because you wouldn't
#   necessarily excise all transcripts in this intergenic region in phase 2's
#   split! (e.g. when an intergenic transcript overlaps with only one adjacent
#   gene's UTR).
# Let's store a list of intergenic regions that are going to be split.
require 'set'
intergenics_to_split = {}
transcripts.transcripts_to_split.each do |chromosome, transcripts_to_split_by_chromosome|
  intergenics_to_split[chromosome] = Set.new
  transcripts_to_split_by_chromosome.each do |_, split_events|
    split_events.each do |split_event|
      intergenics_to_split[chromosome].add split_event
    end
  end
end

# First, let's split transcripts over multiple genes, so that we can assign
#   correct gene_ids to each transcript fragment.
transcripts.split!

# Flag adjacent upstream and downstream transcripts.
puts "#{Time.new}: fixing additional transcripts in these intergenic regions."
intergenics_to_split.each do |chromosome, intergenics_to_split_by_chromosome|
  intergenics_to_split_by_chromosome.each do |split_event|
    transcripts.transcript_ids_for_gene(chromosome, split_event\
        [:upstream_gene_id]).concat(transcripts.transcript_ids_for_gene( \
        chromosome, split_event[:downstream_gene_id])).each do \
        |adjacent_transcript_id|
      if adjacent_transcript_id
        transcripts.define_split(chromosome, adjacent_transcript_id, \
            split_event[:coords], split_event[:upstream_gene_id], \
            nil) # Don't need downstream_gene_id here.
      end
    end
  end
end

transcripts.split!

# Delete intergenic transcripts.
intergenics_to_split.each do |chromosome, intergenics_to_split_by_chromosome|
  intergenics_to_split_by_chromosome.each do |split_event|
    transcripts.add_event(:phase_one_intergenic, [:intergenics, 1])
    transcripts_deleted = transcripts.delete_transcripts_for_gene(chromosome, \
        [split_event[:upstream_gene_id], split_event[:downstream_gene_id]])
    transcripts.add_event(:phase_one_intergenic, [:transcripts, \
        transcripts_deleted])
  end
end

transcripts.split!

# Make (possibly large) object available for garbage collection.
intergenics_to_split = nil

################################################################################
### Find overlapping transcripts with different gene IDs
# TODO: what about ALAD-SPP? Only terminal exons? or not if they share the same tss? Getting complicated!
# If in-gene transcripts overlap with those on adjacent genes, or within
#   adjacent intergenic regions, excise the overlap from the in-gene transcripts
#   and remove the intergenic transcripts.

if !$options[:minimal_split]
  puts "#{Time.new}: re-sorting transcripts."
  transcripts.sort!

  puts "#{Time.new}: fixing overlapping transcripts on adjacent genes."
  refgff.chromosome_names.each do |chromosome|
    if $options[:verbosity] >= 1
      puts "#{Time.new}:   analysing adjacent overlapping transcripts for #{chromosome}."
    end
    prev_gene_transcripts_and_coords = nil
    prev_gene_id = nil
    refgff.gene_ids(chromosome).each do |current_gene_id|
      current_gene_transcripts_and_coords = transcripts.
          transcripts_and_coords_union_for_gene(chromosome, current_gene_id)
      intergenic_transcripts_and_coords = transcripts.
          transcripts_and_coords_union_for_gene(chromosome, \
          [prev_gene_id, current_gene_id])
      if prev_gene_id # Not the first iteration.
          # (Given there are any transcripts for that region,) is there direct
          #   overlap between the in-gene transcripts, or between any in-gene
          #   transcript and an intergenic transcript.
        to_split = false
        if (prev_gene_transcripts_and_coords && \
            current_gene_transcripts_and_coords) && \
            (prev_gene_transcripts_and_coords.last.last >= \
            current_gene_transcripts_and_coords.last.first)
          if $options[:verbosity] >= 2
            puts 'there is overlap between transcript(s) on genes ' \
                  "#{prev_gene_id} and #{current_gene_id}"
          end
          to_split = true
          transcripts.add_event(:phase_two, [:intergenics_adj_genes, 1])
        elsif (prev_gene_transcripts_and_coords && \
            intergenic_transcripts_and_coords) && \
            (prev_gene_transcripts_and_coords.last.last >= \
            intergenic_transcripts_and_coords.last.first)
          if $options[:verbosity] >= 2
            puts 'there is overlap between transcript(s) on gene ' \
                  "#{prev_gene_id} and the downstream intergenic region"
          end
          to_split = true
          transcripts.add_event(:phase_two, [:intergenics_adj_inter, 1])
        elsif (current_gene_transcripts_and_coords && \
            intergenic_transcripts_and_coords) && \
            (intergenic_transcripts_and_coords.last.last >= \
            current_gene_transcripts_and_coords.last.first)
          if $options[:verbosity] >= 2
            puts 'there is overlap between transcript(s) on gene ' \
                  "#{current_gene_id} and the upstream intergenic region"
          end
          to_split = true
          transcripts.add_event(:phase_two, [:intergenics_adj_genes, 1])
        end
        if to_split
          split_coords = [[intergenic_transcripts_and_coords.last.first, \
              current_gene_transcripts_and_coords.last.first].min,
              [prev_gene_transcripts_and_coords.last.last, \
              intergenic_transcripts_and_coords.last.last].max]
          # Cut prev_gene and current_gene transcripts.
          (prev_gene_transcripts_and_coords.first.keys + \
              current_gene_transcripts_and_coords.first.keys).each do |transcript_id|
            transcripts.define_split(chromosome, transcript_id, split_coords, \
                  prev_gene_id, nil) # Don't need downstream_gene_id here.
            transcripts.write_gene_id(chromosome, transcript_id, current_gene_id)
          end
          # Delete intergenic transcripts.
          transcripts.add_event(:phase_two, [:transcripts, \
              transcripts.delete_transcripts_for_gene(chromosome, \
              [prev_gene_id, current_gene_id])])
        end
      end
      prev_gene_transcripts_and_coords = current_gene_transcripts_and_coords
      prev_gene_id = current_gene_id
    end
  end

  transcripts.split!
end

################################################################################
### Unimplemented: Make gene name unique if transcripts do not overlap
# DEXSeq trusts geneIDs. Hence, it combines two genes if they have the same
#   geneID, regardless of where they are located.
# TODO: intergenic transcripts should have unique gene names: e.g. *:a, *:b.
#   Not sure if I'll bother with this. Our premise is to trust the gene models.
#   Allowing intergenic transcripts to exist at all is beyond our premise.
# N.B. transcripts within a single gene that don't overlap will still have the
#   same gene ID. I won't change this. Again, trust the gff.
################################################################################

transcripts.write_to_file($options[:output_path])

puts "#{Time.new}: processing complete and file written.

SUMMARY
=======
Analysis of transcripts that overlap multiple genes:
"

phase_one_overlaps_stats = transcripts.stats[:phase_one_overlaps]
total_transcripts = 0
phase_one_overlaps_stats.values.each do |count|
  total_transcripts += count
end

if (na_count = phase_one_overlaps_stats[:NA])
  puts "  #{na_count} transcript#{'s' if na_count != 1} (" \
      "#{(na_count.to_f/total_transcripts*100).round(1)}%) "\
      'with no matching reference contig'
  phase_one_overlaps_stats.delete(:NA)
end
phase_one_overlaps_stats.sort_by {|key, _| key }.each do |genes, count|
  puts "  #{count} transcript#{'s' if count != 1} (" \
      "#{(count.to_f/total_transcripts*100).round(1)}%) overlap with #{genes} " \
      "gene#{'s' if genes != 1} each"
end

puts "  ------------------------
  #{total_transcripts} transcripts total

  #{transcripts.stats[:phase_one_intergenic][:intergenics]} intergenic regions " \
"removed, including #{transcripts.stats[:phase_one_intergenic][:transcripts]} " \
"associated intergenic transcripts.\n\n"

if $options[:minimal_split]
  puts 'Minimal-split option specified. Overlapping transcripts on adjacent ' \
      'genes have been ignored.'
else
  puts "Transcripts on adjacent genes that overlap:
  #{transcripts.stats[:phase_two][:intergenics_adj_genes] + \
  transcripts.stats[:phase_two][:intergenics_adj_inter]} intergenic regions " \
  "fixed.
    #{transcripts.stats[:phase_two][:intergenics_adj_genes]} on " \
    "adjacent genes (possibly including intergenics).
    #{transcripts.stats[:phase_two][:intergenics_adj_inter]} purely from " \
    "overlaps with intergenic transcripts.
  #{transcripts.stats[:phase_two][:transcripts]} associated intergenic " \
  'transcripts removed.'
end