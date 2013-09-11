all_events = new.env()
for (j in 1:5) {
	all_events_method = new.env()

	for (index in 1:5) {
		if (index == j) 
			next
		method = methods[index]
		data = read.table(sprintf("%s/final/%s", inputDir, method))
		chrs = unique(data[, 1])
		for (chr in chrs) {
			nextInterval = (Intervals(data[which(data[, 1] == chr), 2:3]))
			all_events_method[[chr]] = (Intervals(rbind(all_events_method[[chr]], nextInterval)))
		}
	}
	all_events[[sprintf("%d", j)]] = all_events_method
}

getNumOverlaps = function(int, all) {
	if (is.null(all)) 
		return(0)
	all = reduce(all)
	count = 0
	if (dim(int)[1] == 0) 
		return(0)
	overlaps = interval_overlap(int, all)
	sumoverlap = 0
	for (i in unlist(overlaps)) {
		candidate = all[i]
		sumoverlap = sumoverlap + min(candidate[, 2], int[, 2]) - max(candidate[, 1], int[, 1])
	}
	if (sumoverlap > 0.5 * size(int)) 
		count = count + 1
	return(count)
}

##### This part shows the plot of number of events who have a set of corresponding other method events that cover at least 50% of this event
for (index in 1:5) {
	method = methods[index]
	label = labels[index]
	color = colors[index]
	data = read.table(sprintf("%s/final/%s", inputDir, method))
	data = data[which(data[, 4] == "no"), ]
	event = new.env()
	chrs = unique(data[, 1])

	for (chr in chrs) event[[chr]] = Intervals(data[which(data[, 1] == chr), 2:3])

	starts = c(0, 50, 200, 400, 600, 800, 1000, 3000, 5000, 7000, 9000, 20 * 10^3, 200 * 10^3, 2000 * 10^3, 20000 * 10^3, 2e+05 * 10^3)
	ends = c(50, 200, 400, 600, 800, 1000, 3000, 5000, 7000, 9000, 20 * 10^3, 200 * 10^3, 2000 * 10^3, 20000 * 10^3, 2e+05 * 10^3, 2 * 109)
	count = rep(0, length(starts))
	count2 = count
	count3 = count
	for (chr in chrs) {
		evs = event[[chr]]
		if (is.null(evs)) 
			next
		for (j in 1:dim(evs)[1]) {
			ev = evs[j, ]
			if (is.null(ev)) 
				next
			for (i in 1:length(starts)) {
				if (size(ev) >= starts[i] && size(ev) < ends[i]) {
					over = getNumOverlaps(ev, all_events[[sprintf("%d", index)]][[chr]])
					count[i] = count[i] + (1 + over)
					if (over > 1) 
						over = 1
					count3[i] = count3[i] + (1 + over)
					count2[i] = count2[i] + 1
				}
			}
		}
	}
	png(file = sprintf("%sdouble_multi_no_plots/%s.png", inputDir, label), width = 900, height = 900)
	par(ps = 24, cex = 1, cex.main = 1, cex.lab = 1, cex.sub = 1, cex.axis = 1)
	colors = c("red", "forestgreen", "black", "darkblue", "cyan", "yellow")
	plot((starts + ends)/2, count, type = "b", log = "x", col = color, lwd = 6, lty = 1, xaxt = "n")
	lines((starts + ends)/2, count2, col = "purple", lwd = 6, lty = 1, type = "b")
	lines((starts + ends)/2, count3, col = "maroon", lwd = 6, lty = 1, type = "b")

	axis(1, at = (starts + ends)/2, labels = (starts + ends)/2)
	legend("topleft", col = color, legend = label, lwd = 2, lty = 1)
	dev.off()
}