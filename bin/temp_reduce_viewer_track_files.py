infile = '/Users/michaelstadler/Bioinformatics/Projects/insulators/viewerfiles/viewer_tracks/synteny-breaks-m-v.txt'
outfile = '/Users/michaelstadler/Bioinformatics/Projects/insulators/viewerfiles/viewer_tracks_reduced/synteny-breaks-m-v.txt'

out = open(outfile, 'w')
in_ = open(infile, 'r')
for line in in_:
	line = line.rstrip()
	chr_, pos1, pos2, val = line.split()
	bin_ = int(int(pos1) / 500)
	val = float(val)
	if val != 0:
		out.write(chr_ + '\t' + str(bin_) + '\t' + "{:.2f}".format(val) + '\n')

out.close()
in_.close()