###### This is a simple table that shows what fraction of the events of one method overlap with other events of the same method
table = NULL
for (index in 1:5) {
	method = methods[index]
	label = labels[index]
	color = colors[index]
	data = read.table(sprintf("%s/final/%s", inputDir, method))
	#data=data[which(data[,4]=="no"),]
	event = new.env()
	reduced_event = new.env()
	chrs = unique(data[, 1])
	sum_e = 0
	sum_re = 0
	for (chr in chrs) {
		event[[chr]] = Intervals(data[which(data[, 1] == chr), 2:3])
		reduced_event[[chr]] = reduce(event[[chr]])
		sum_e = sum_e + dim(event[[chr]])[1]
		sum_re = sum_re + dim(reduced_event[[chr]])[1]
	}
	table = rbind(table, cbind(method, sum_e, sum_re, sum_re/sum_e))
}

###### This is a similar table which shows which events overlap with the segmental duplication regions
seg_dups = read.table(sprintf("%s/seg_dups.bed", inputDir))
seg_dup_int = new.env()
for (chr in chrs) seg_dup_int[[chr]] = Intervals(seg_dups[which(seg_dups[, 1] == chr), 2:3])

table = NULL
for (index in 1:5) {
	method = methods[index]
	label = labels[index]
	color = colors[index]
	data = read.table(sprintf("%s/final/%s", inputDir, method))
	data = data[which(data[, 4] == "asm"), ]
	event = new.env()
	reduced_event = new.env()
	chrs = unique(data[, 1])
	sum_e = 0
	sum_re = 0
	for (chr in chrs) {
		event[[chr]] = Intervals(data[which(data[, 1] == chr), 2:3])
		if (is.null(seg_dup_int[[chr]]) == FALSE) {
			for (t in 1:dim(event[[chr]])[1]) {
				next_ev = event[[chr]][i, ]
				reduced_event[[chr]] = interval_difference(next_ev, seg_dup_int[[chr]])
				sum_e = sum_e + 1
				if (dim(reduced_event[[chr]])[1] > 0) 
					sum_re = sum_re + 1
			}
		} else {
			reduced_event[[chr]] = event[[chr]]
		}
	}
	table = rbind(table, cbind(method, sum_e, sum_re, sum_re/sum_e))
}
