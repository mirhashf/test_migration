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
		if (min(candidate[, 2], int[, 2]) - max(candidate[, 1], int[, 1]) > 0.5 * size(int)) 
			#if(size(interval_intersection(int,candidate)>0.5*size(int)))
			count = count + 1
	}
	return(count)
}



##### Generates the mappability score plots for all assembled events
for (index in 1:5) {
	colors = c("red", "forestgreen", "black", "darkblue", "cyan", "yellow")
	inputDir = "~/expvalid_size/"
	methods = dir(sprintf("%s/final/", inputDir))
	method = methods[index]
	maps = dir(sprintf("%s/map/", inputDir))
	map = maps[index]
	labels = c("1KG", "Complete", "Conrad", "Kidd", "McCaroll")
	label = labels[index]
	color = colors[index]
	data = read.table(sprintf("%s/final/%s", inputDir, method))
	data = data[which(data[, 4] == "asm"), ]
	map_data = read.table(sprintf("%s/map/%s", inputDir, map))
	event = new.env()
	map_event = new.env()
	chrs = unique(data[, 1])
	for (chr in chrs) {
		event[[chr]] = Intervals(data[which(data[, 1] == chr), 2:3])
		u = map_data[which(data[, 1] == chr), 4:5]
		u[(is.na(u))] <- 0
		map_event[[chr]] = u
	}
	starts = c(0, 50, 200, 400, 600, 800, 1000, 3000, 5000, 7000, 9000, 20 * 10^3, 200 * 10^3, 2000 * 10^3, 20000 * 10^3, 2e+05 * 10^3)
	ends = c(50, 200, 400, 600, 800, 1000, 3000, 5000, 7000, 9000, 20 * 10^3, 200 * 10^3, 2000 * 10^3, 20000 * 10^3, 2e+05 * 10^3, 2 * 109)
	count0 = rep(0, length(starts))
	count = rep(0, length(starts))
	count2 = count
	count3 = count
	for (chr in chrs) {
		evs = event[[chr]]
		for (j in 1:dim(evs)[1]) {
			ev = evs[j, ]
			for (i in 1:length(starts)) {
				if (size(ev) >= starts[i] && size(ev) < ends[i]) {
					range = (map_event[[chr]])[j, ]
					lower = range[1] - range[2]
					upper = range[1] + range[2]
					count0[i] = count0[i] + as.numeric(range[1])
					count[i] = count[i] + as.numeric(lower)
					count3[i] = count3[i] + as.numeric(upper)
					count2[i] = count2[i] + 1
				}
			}
		}
	}
	count = count/count2
	count[which(is.nan(count))] = 0
	count3 = count3/count2
	count3[which(is.nan(count3))] = 0
	count0 = count0/count2
	count0[which(is.nan(count0))] = 0
	png(file = sprintf("%s/mappability_asm/%s.png", inputDir, label), width = 900, height = 900)
	par(ps = 24, cex = 1, cex.main = 1, cex.lab = 1, cex.sub = 1, cex.axis = 1)
	colors = c("red", "forestgreen", "black", "darkblue", "cyan", "yellow")
	plot((starts + ends)/2, count0, type = "b", log = "x", col = color, lwd = 6, lty = 1, xaxt = "n")
	lines((starts + ends)/2, count3, type = "b", log = "x", col = color, lwd = 3, lty = 3)
	lines((starts + ends)/2, count, type = "b", log = "x", col = color, lwd = 3, lty = 3)
	axis(1, at = (starts + ends)/2, labels = (starts + ends)/2)
	legend("topleft", col = color, legend = label, lwd = 2, lty = 1)
	dev.off()
}



### This one generates the above plot for only those events whose mappability score is greater than the threshold
thresh = 0.75
for (index in 1:5) {
	colors = c("red", "forestgreen", "black", "darkblue", "cyan", "yellow")
	inputDir = "~/expvalid_size/"
	methods = dir(sprintf("%s/final/", inputDir))
	method = methods[index]
	maps = dir(sprintf("%s/map/", inputDir))
	map = maps[index]

	labels = c("1KG", "Complete", "Conrad", "Kidd", "McCaroll")
	label = labels[index]
	color = colors[index]
	map_data = read.table(sprintf("%s/map/%s", inputDir, map))
	indices = which(map_data[, 4] > thresh)
	data = read.table(sprintf("%s/final/%s", inputDir, method))
	data = data[indices, ]
	data = data[which(data[, 4] == "asm"), ]

	event = new.env()
	map_event = new.env()
	chrs = unique(data[, 1])

	for (chr in chrs) {
		event[[chr]] = Intervals(data[which(data[, 1] == chr), 2:3])
	}
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
	filename = sprintf("%s/mappability_double_plots/%g/", inputDir, thresh)
	dir.create(filename, recursive = TRUE)
	png(file = sprintf("%s/%s.png", filename, label), width = 900, height = 900)
	par(ps = 24, cex = 1, cex.main = 1, cex.lab = 1, cex.sub = 1, cex.axis = 1)
	colors = c("red", "forestgreen", "black", "darkblue", "cyan", "yellow")
	plot((starts + ends)/2, count, type = "b", log = "x", col = color, lwd = 6, lty = 1, xaxt = "n")
	lines((starts + ends)/2, count2, col = "pink", lwd = 6, lty = 2, type = "b")
	lines((starts + ends)/2, count3, col = "yellow", lwd = 6, lty = 2, type = "b")

	axis(1, at = (starts + ends)/2, labels = (starts + ends)/2)
	legend("topleft", col = color, legend = label, lwd = 2, lty = 1)
	dev.off()

}

##### This is an important function that generates the "validated" table for all events
all_event = NULL
validates = NULL
num_over = NULL
for (index in 1:5) {
	detail_event = NULL
	inputDir = "~/expvalid_size/"
	methods = dir(sprintf("%s/final/", inputDir))
	method = methods[index]
	maps = dir(sprintf("%s/map/", inputDir))
	map = maps[index]
	labels = c("1KG", "Complete", "Conrad", "Kidd", "McCaroll")
	label = labels[index]
	color = colors[index]
	data = read.table(sprintf("%s/final/%s", inputDir, method))
	map_data = read.table(sprintf("%s/map/%s", inputDir, map))
	event = new.env()
	map_event = new.env()
	chrs = unique(data[, 1])
	for (chr in chrs) {
		u = map_data[which(data[, 1] == chr), 4:5]
		u[(is.na(u))] <- 0
		map_event[[chr]] = u
		event[[chr]] = Intervals(data[which(data[, 1] == chr), 2:3])
		detail_event = rbind(detail_event, cbind(label, chr, data[which(data[, 1] == chr), 2:3], map_event[[chr]], data[which(data[, 1] == chr), 4]))
	}
	validates = rep("no", dim(detail_event)[1])
	num_over = rep(1, dim(detail_event)[1])
	u = which(detail_event[, 1] == "Complete" & detail_event[, 7] == "asm" & detail_event[, 5] > 0.75)
	validates[u] = "yes"
	u = which(detail_event[, 7] == "yes")
	validates[u] = "yes"

	starts = c(0, 50, 200, 400, 600, 800, 1000, 3000, 5000, 7000, 9000, 20 * 10^3, 200 * 10^3, 2000 * 10^3, 20000 * 10^3, 2e+05 * 10^3)
	ends = c(50, 200, 400, 600, 800, 1000, 3000, 5000, 7000, 9000, 20 * 10^3, 200 * 10^3, 2000 * 10^3, 20000 * 10^3, 2e+05 * 10^3, 2 * 109)
	count0 = rep(0, length(starts))
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
					if (over > 1) {
						ind = which(detail_event[, 1] == label & detail_event[, 2] == chr & detail_event[, 3] == ev[, 1] & detail_event[, 4] == ev[, 2])
						if (length(ind) != 1) 
							print(ind)
						validates[ind] = "yes"
						num_over[ind] = over
					}
				}
			}
		}
	}
	detail_event = cbind(detail_event, validates, num_over)
	all_event = rbind(all_event, detail_event)
}
names(all_event) = c("Method", "Chr", "Start", "End", "MeanMapScore", "SdMapScore", "ExpValid", "Validated", "NumOverlap")
write.table(all_event, file = sprintf("%s/final_report.txt", inputDir), quote = FALSE, row.names = FALSE, col.names = TRUE)
