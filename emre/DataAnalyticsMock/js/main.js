$(document).ready(function() {
  $("#header").html('<div class="navbar navbar-default navbar-fixed-top">' + 
		'<div class="container">' +
		'<div class="navbar-header">' +
			'<button type="button" class="navbar-toggle" data-toggle="collapse" data-target=".navbar-collapse">' +
				'<span class="icon-bar"></span> <span class="icon-bar"></span> <span class="icon-bar"></span>' +
			'</button>' +
			'<a class="navbar-brand" href="#"><img src="image/leaf-brand.png" style="margin-right: 10px;">Bina Analytics</a>' +
		'</div>' +
		'<div class="collapse navbar-collapse">' +
			'<ul class="nav navbar-nav">' +
				'<li id="data_tab"><a href="data.html"><span class="glyphicon glyphicon-hdd"></span> Data</a></li>' +
				'<li id="datasets_tab"><a href="datasets.html"><span class="glyphicon glyphicon-list-alt"></span> Data sets</a></li>' +
				'<li id="analyses_tab" class="dropdown"><a href="" class="dropdown-toggle" data-toggle="dropdown"><span class="glyphicon glyphicon-random"></span> Analyses <b class="caret"></b></a>' +
					'<ul class="dropdown-menu">' +
						'<li><a href="annotation_analysis_select_sample.html">Annotation</a></li>' +
						'<li><a href="trio_analysis_select_trio.html">Trio Analysis</a></li>' +
						'<li class="disabled"><a href="case_control_study.html">Case/Control study</a></li>' +
					'</ul></li>' +
				'<li id="saved_queries_tab"><a href="saved_queries.html"><span class="glyphicon glyphicon-pushpin"></span> Saved Queries</a></li>' +
				'<li id="ann_sources_tab"><a href="annotation_sources.html"><span class="glyphicon glyphicon-folder-open"></span> Ann. sources</a></li>' +
			'</ul>' +
			'<ul class="nav navbar-nav navbar-right">' +
		      '<li class="dropdown">' +
		        '<a href="#" class="dropdown-toggle" data-toggle="dropdown"><span class="glyphicon glyphicon-user"></span> Hello Dude! <b class="caret"></b></a>' +
		        '<ul class="dropdown-menu">' +
		          '<li><a href="#">Settings</a></li>' +
		          '<li><a href="#">Logout</a></li>' +
		        '</ul>' +
		      '</li>' +
		    '</ul>' +
		'</div>' +
	'</div>' +
	'</div>');
  
  $("#footer").html('<div class="container"><p class="text-muted credit"><a href="http://binatechnologies.com">Bina Technologies Inc.</a></p></div>');
});