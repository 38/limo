#!/uufs/chpc.utah.edu/common/home/u0875014/.linuxbrew/Cellar/python/3.6.5_1/bin/python3.6
import sys
size = {}
size["1"] = 249250621
size["2"] = 243199373
size["3"] = 198022430
size["4"] = 191154276
size["5"] = 180915260
size["6"] = 171115067
size["7"] = 159138663
size["8"] = 146364022
size["9"] = 141213431
size["10"] = 135534747
size["11"] = 135006516
size["12"] = 133851895
size["13"] = 115169878
size["14"] = 107349540
size["15"] = 102531392
size["16"] = 90354753
size["17"] = 81195210
size["18"] = 78077248
size["19"] = 59128983
size["20"] = 63025520
size["21"] = 48129895
size["22"] = 51304566
size["X"] = 155270560
size["Y"] = 59373566

last_chrom = "1"
last_end = 0

for line in sys.stdin:
	line = line.strip().split()
	if line[0] not in size: 
		continue
	line[1] = int(line[1])
	line[2] = int(line[2])
	if last_chrom != line[0]:
		print(last_chrom, last_end, size[last_chrom], sep = "\t")
		last_end = 0
		last_chrom = line[0]

	if last_end < line[1]:
		print(line[0], last_end, line[1], sep="\t")
	last_end = max(last_end, line[2])

if last_end != size[last_chrom]:
	print(last_chrom, last_end, size[last_chrom], sep = "\t")
	
