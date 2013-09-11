#### This part generates the number events of particular size for each method
for (index in 1:5) {
	colors = c("red", "forestgreen", "black", "darkblue", "cyan", "yellow")

	inputDir = "~/expvalid_size/"
	methods = dir(sprintf("%s/final/", inputDir))
	method = methods[index]
	labels = c("1KG", "Complete", "Conrad", "Kidd", "McCaroll")
	label = labels[index]
	color = colors[index]
	data = read.table(sprintf("%s/final/%s", inputDir, method))
	event = new.env()
	chrs = unique(data[, 1])

	for (chr in chrs) event[[chr]] = Intervals(data[which(data[, 1] == chr), 2:3])

	starts = c(0, 50, 200, 400, 600, 800, 1000, 3000, 5000, 7000, 9000, 20 * 10^3, 200 * 10^3, 2000 * 10^3, 20000 * 10^3, 2e+05 * 10^3)
	ends = c(50, 200, 400, 600, 800, 1000, 3000, 5000, 7000, 9000, 20 * 10^3, 200 * 10^3, 2000 * 10^3, 20000 * 10^3, 2e+05 * 10^3, 2 * 109)
	count = rep(0, length(starts))

	for (chr in chrs) {
		evs = event[[chr]]
		for (j in 1:dim(evs)[1]) {
			ev = evs[j, ]
			for (i in 1:length(starts)) {
				if (size(ev) >= starts[i] && size(ev) < ends[i]) 
					count[i] = count[i] + 1
			}
		}
	}
	png(file = sprintf("%s/single_plots/%s.png", inputDir, label), width = 900, height = 900)
	par(ps = 24, cex = 1, cex.main = 1, cex.lab = 1, cex.sub = 1, cex.axis = 1)
	colors = c("red", "forestgreen", "black", "darkblue", "cyan", "yellow")
	plot((starts + ends)/2, count, type = "b", log = "x", col = color, lwd = 6, lty = 1, xaxt = "n")
	axis(1, at = (starts + ends)/2, labels = (starts + ends)/2)
	legend("topleft", col = color, legend = label, lwd = 2, lty = 1)
	dev.off()
}

all_events = new.env()

for (index in 1:5) {
	method = methods[index]
	data = read.table(sprintf("%s/final/%s", inputDir, method))
	chrs = unique(data[, 1])
	for (chr in chrs) {
		nextInterval = Intervals(data[which(data[, 1] == chr), 2:3])
		all_events[[chr]] = Intervals(rbind(all_events[[chr]], nextInterval))
	}
}

getNumOverlaps = function(int, all) {
	count = 0
	overlaps = interval_overlap(int, all)
	for (i in unlist(overlaps)) {
		candidate = all[i]
		#if(min(candidate[,2],int[,2])-max(candidate[,1],int[,1])>0.999*size(int))
		if (size(interval_intersection(int, candidate) > 0.5 * size(int))) 
			count = count + 1
	}
	return(count)
}
##### This part shows for each event, how many events of other methods overlap with it.

for (index in 1:5) {
	method = methods[index]
	label = labels[index]
	color = colors[index]
	data = read.table(sprintf("%s/final/%s", inputDir, method))
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
		for (j in 1:dim(evs)[1]) {
			ev = evs[j, ]
			for (i in 1:length(starts)) {
				if (size(ev) >= starts[i] && size(ev) < ends[i]) {
					over = getNumOverlaps(ev, all_events[[chr]])
					count[i] = count[i] + over
					if (over > 2) 
						over = 2
					count3[i] = count3[i] + over
					count2[i] = count2[i] + 1
				}
			}
		}
	}
	png(file = sprintf("%sdouble_plots/%s.png", inputDir, label), width = 900, height = 900)
	par(ps = 24, cex = 1, cex.main = 1, cex.lab = 1, cex.sub = 1, cex.axis = 1)
	colors = c("red", "forestgreen", "black", "darkblue", "cyan", "yellow")
	plot((starts + ends)/2, count, type = "b", log = "x", col = color, lwd = 6, lty = 1, xaxt = "n")
	lines((starts + ends)/2, count2, col = "pink", lwd = 6, lty = 2, type = "b")
	lines((starts + ends)/2, count3, col = "yellow", lwd = 6, lty = 2, type = "b")
	axis(1, at = (starts + ends)/2, labels = (starts + ends)/2)
	legend("topleft", col = color, legend = label, lwd = 2, lty = 1)
	dev.off()
}
