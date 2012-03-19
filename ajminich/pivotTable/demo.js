
var header = ["SeqAlto Alignment","BWA Alignment","Quantity","Value"];
var options = ["Mapped","Multiple-Mapped","Unmapped"];
var quantities = ["Number of Reads","Gap Opens","Mismatches"]


// Create data
var data = [];
for (var i=0; i < options.length; i++)
for (var j=0; j < options.length; j++)
for (var k=0; k < quantities.length; k++) {
	var value = Math.round(Math.random()*100);
	data.push([options[i],options[j],quantities[k],value]);
}

var pivot = {
	panel:1,
	tab:5,
	div:"pivot",
	needs:["pivot","statistics"],
	cb:function() {
		var pivot = new OAT.Pivot("pivot_content","pivot_chart","pivot_page",header,data,[0],[1],[2],3,{showChart:1});
		
		var aggRef = function() {
			pivot.options.agg = parseInt($v("pivot_agg"));
			pivot.go();
		}
		var aggRefTotals = function() {
			pivot.options.aggTotals = parseInt($v("pivot_agg_totals"));
			pivot.go();
		}
		
		/* create agg function list */
		OAT.Dom.clear("pivot_agg");
		OAT.Dom.clear("pivot_agg_totals");
		for (var i=0;i<OAT.Statistics.list.length;i++) {
			var item = OAT.Statistics.list[i];
			OAT.Dom.option(item.shortDesc,i,"pivot_agg");
			OAT.Dom.option(item.shortDesc,i,"pivot_agg_totals");
			if (pivot.options.agg == i) { $("pivot_agg").selectedIndex = i; }
			if (pivot.options.aggTotals == i) { $("pivot_agg_totals").selectedIndex = i; }
		}
		OAT.Dom.attach("pivot_agg","change",aggRef);
		OAT.Dom.attach("pivot_agg_totals","change",aggRefTotals);
	}
}

function init() {

	if (pivot.cb) {
		var ref = function() {
			pivot.cb();
			pivot.drawn = true;
		}
		OAT.Loader.loadFeatures(pivot.needs,ref);
	}

}
