# Convert gff3 to BED12 format
# Usage: julia gff3ToBed12.jl path_gff3 path_bed


# Assumption:
# - The input gff3 file has to contain rows of "exon", "start_codon", "stop_codon".
# - chromStrat and chromEnd is set to the left-most exon's start (0-based) and the right-most exon's end positions, respectively.
# - thickStart and thickEnd is set to the minimum (0-based) and maximum (1-based) coordinates of start_codon and stop_codon.
# - If a transcript_id doesn't have start_codon or stop codon, thickStart and thickEnd is set to the value of chromEnd.

# Test:
	# $ grep ENSMUST00000070533 Mouse_UCSC_mm10.Ensembl.bed
# 	# chr1    3214481 3671498 ENSMUST00000070533      0       -       3216021 3671348 0       3       2487,200,947,   0,207220,456070,
	# $ grep ENSMUST00000070533 gencode.vM9.annotation.bed
	# chr1    3214481 3671498 ENSMUST00000070533.4    0       -       3216021 3671348 0       3       2487,200,947,   0,207220,456070,

function gff3ToBed12(path_gff3, path_bed)
	transcript_id = ""
	chr = ""
	set_start = BitSet()
	set_end = BitSet()
	set_codon = BitSet()
	strand = ""

	fw = open(path_bed, "w")

	f = open(path_gff3)
	for ln in eachline(f)
		if occursin(r"^#", ln)
			continue
		end

		arr = split(chomp(ln), "\t")
		if arr[3] == "gene"
			continue
		end

		new_transcript_id = match(r"transcript_id=(.+?);", arr[9]).captures[1]
		if transcript_id != new_transcript_id
			if length(set_start) > 0
				res = calcLine(chr, set_start, set_end, transcript_id, strand, set_codon)
				write(fw, res * "\n")
			end

			transcript_id = new_transcript_id
			chr = arr[1]
			empty!(set_start)
			empty!(set_end)
			empty!(set_codon)
			strand = arr[7]
		end

		if arr[3] == "exon"
			push!(set_start, parse(Int64, arr[4]))
			push!(set_end, parse(Int64, arr[5]))
		elseif arr[3] == "stop_codon" || arr[3] == "start_codon"
			push!(set_codon, parse(Int64, arr[4]))
			push!(set_codon, parse(Int64, arr[5]))
		end

	end
	close(f)


	if length(set_start) > 0
		res = calcLine(chr, set_start, set_end, transcript_id, strand, set_codon)
		write(fw, res * "\n")
	end
	close(fw)
end

function calcBlockSizes(set_start::BitSet, set_end::BitSet)
	join([string(collect(set_end)[i] - collect(set_start)[i] + 1) * "," for i in range(1, stop=length(set_start))], "")
end

function calcBlockStarts(set_start::BitSet, chromStart)
	join([string(i - chromStart) * "," for i in set_start ], "")
end

function calcThickStart(set_codon::BitSet, chromEnd)
	if length(set_codon) == 0
		chromEnd
	else
		first(set_codon) - 1
	end
end

function calcThickEnd(set_codon::BitSet, chromEnd)
	if length(set_codon) == 0
		chromEnd
	else
		last(set_codon)
	end
end

function calcLine(
					chr::AbstractString, set_start::BitSet, set_end::BitSet,
					transcript_id::AbstractString, strand::AbstractString,
					set_codon::BitSet
					)
	chrom = chr
	if length(set_start) == 0
		print(transcript_id)
	end

	chromStart = first(set_start) - 1
	chromEnd = last(set_end)
	name = transcript_id
	score = 0
	strand = strand
	thickStart = calcThickStart(set_codon, chromEnd)
	thickEnd = calcThickEnd(set_codon, chromEnd)
	itemRgb = 0
	blockCount = length(set_start)
	blockSizes = calcBlockSizes(set_start, set_end)
	blockStarts = calcBlockStarts(set_start, chromStart + 1)

	res = join([
				chrom;
				string(chromStart);
				string(chromEnd);
				name;
				string(score);
				strand;
				string(thickStart);
				string(thickEnd);
				string(itemRgb);
				blockCount;
				blockSizes;
				blockStarts;
			], "\t")
	return res
end

# chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
# chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
# chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
# name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
# score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray: shade score in range  	≤ 166	167-277	278-388	389-499	500-611	612-722	723-833	834-944	≥ 945
# strand - Defines the strand - either '+' or '-'.
# thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
# thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
# itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
# blockCount - The number of blocks (exons) in the BED line.
# blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
# blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.

path_gff3 = ARGS[1]
path_bed = ARGS[2]


gff3ToBed12(path_gff3, path_bed)
