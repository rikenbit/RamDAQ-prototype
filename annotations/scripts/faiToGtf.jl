function faiToGtf(pathFai, pathGtf; addLen = 0)
    mkpath(dirname(pathGtf))

    fw = open(pathGtf, "w")
    open(pathFai) do f
        for ln in eachline(f)
            arr = split(chomp(ln), "\t")
            oldStart = 0
            oldEnd = parse(Int64, arr[2])
            newStart = string(oldStart + 1 - addLen)
            newEnd =   string(oldEnd + addLen)
            write(fw, @sprintf "%s\tbedToGtf\texon\t%s\t%s\t%s\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"\n" arr[1] newStart newEnd "0" "." arr[1] arr[1])

        end
    end
    close(fw)

    # 1. seqname - The name of the sequence. Must be a chromosome or scaffold.
    # 2. source - The program that generated this feature.
    # 3. feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
    # 4. start - The starting position of the feature in the sequence. The first base is numbered 1.
    # 5. end - The ending position of the feature (inclusive).
    # 6. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
    # 7 s trand - Valid entries include '+', '-', or '.' (for don't know/don't care).
    # 8 frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
    # group - All lines with the same group are linked together into a single item.
end


pathFai = ARGS[1]
pathGtf = ARGS[2]

faiToGtf(pathFai, pathGtf)
