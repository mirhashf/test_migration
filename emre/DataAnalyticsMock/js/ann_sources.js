var ann_source_accordion = 
'<div class="panel-group" id="accordion">' +
	'<div class="panel panel-default">' +
		'<div class="panel-heading">' +
			'<h3 class="panel-title"><a class="accordion-toggle" data-toggle="collapse" data-parent="#accordion" href="#collapse1">Small Variations</a></h3>' +
		'</div>' +
		'<div id="collapse1" class="panel-collapse collapse in">' +
			'<div class="panel-body">' +
				'<div class="row">' +
					'<div class="col-md-3">' +
						'<div class="panel panel-default">' +
							'<div class="panel-heading">' +
								'<h3 class="panel-title"><label> <input type="checkbox" checked> dbSNP (135)</label></h3>' +
							'</div>' +
							'<div class="panel-body ann-source-box">' +
								'<label> <input type="checkbox" checked> name</label> <br/>' +
								'<label> <input type="checkbox" checked> chrom</label> <br/>' +
								'<label> <input type="checkbox" checked> chromStart</label><br/>' +
								'<label> <input type="checkbox" checked> chromEnd</label><br/>' +
								'<label> <input type="checkbox" checked> func</label><br/>' +
								'<label> <input type="checkbox" checked> class</label><br/>' +
								'<label> <input type="checkbox" checked> score</label>' +
							'</div>' +
						'</div>' +
					'</div>' +
		
					'<div class="col-md-3">' +
						'<div class="panel panel-default">' +
							'<div class="panel-heading">' +
								'<h3 class="panel-title"><label> <input type="checkbox" checked> HapMap SNPs</label></h3>' +
							'</div>' +
							'<div class="panel-body ann-source-box">' +
								'<label> <input type="checkbox" checked> chrom</label> <br/>' +
								'<label> <input type="checkbox" checked> chromStart</label><br/>' +
								'<label> <input type="checkbox" checked> chromEnd</label><br/>' +
								'<label> <input type="checkbox" checked> func</label><br/>' +
								'<label> <input type="checkbox" checked> score</label><br/>' +
								'<label> <input type="checkbox" checked> strand</label>' +
							'</div>' +
						'</div>' +
					'</div>' +
					'<div class="col-md-3">' +
						'<div class="panel panel-default">' +
							'<div class="panel-heading">' +
								'<h3 class="panel-title"><label> <input type="checkbox"> 1000 Genomes</label></h3>' +
							'</div>' +
							'<div class="panel-body ann-source-box">' +
								'<label> <input type="checkbox"> col1</label> <br/>' +
								'<label> <input type="checkbox"> col2</label><br/>' +
								'<label> <input type="checkbox"> col3</label><br/>' +
								'<label> <input type="checkbox"> col4</label><br/>' +
								'<label> <input type="checkbox"> col5</label><br/>' +
								'<label> <input type="checkbox"> col6</label>' +
							'</div>' +
						'</div>' +
					'</div>' +
				'</div>' +
			'</div>' +
		'</div>' +
	'</div>' +
	'<div class="panel panel-default">' +
		'<div class="panel-heading">' +
			'<h3 class="panel-title"><a class="accordion-toggle" data-toggle="collapse" data-parent="#accordion" href="#collapse2">Structural variations</a></h3>' +
		'</div>' +
		'<div id="collapse2" class="panel-collapse collapse">' +
			'<div class="panel-body">' +
				'<div class="row">' +
					'<div class="col-md-3">' +
						'<div class="panel panel-default">' +
							'<div class="panel-heading">' +
								'<h3 class="panel-title"><label> <input type="checkbox"> dbVar</label></h3>' +
							'</div>' +
							'<div class="panel-body ann-source-box">' +
								'<label> <input type="checkbox"> col1</label> <br/>' +
								'<label> <input type="checkbox"> col2</label><br/>' +
								'<label> <input type="checkbox"> col3</label><br/>' +
								'<label> <input type="checkbox"> col4</label><br/>' +
								'<label> <input type="checkbox"> col5</label><br/>' +
								'<label> <input type="checkbox"> col6</label>' +
							'</div>' +
						'</div>' +
					'</div>' +
					'<div class="col-md-3">' +
						'<div class="panel panel-default">' +
							'<div class="panel-heading">' +
								'<h3 class="panel-title"><label> <input type="checkbox" checked> DGV</label></h3>' +
							'</div>' +
							'<div class="panel-body ann-source-box">' +
								'<label> <input type="checkbox" checked> chrom</label> <br/>' +
								'<label> <input type="checkbox" checked> chromStart</label><br/>' +
								'<label> <input type="checkbox" checked> chromEnd</label><br/>' +
								'<label> <input type="checkbox" checked> func</label><br/>' +
								'<label> <input type="checkbox" checked> score</label><br/>' +
								'<label> <input type="checkbox" checked> strand</label><br/>' +
								'<label> <input type="checkbox" checked> thickStart</label><br/>' +
								'<label> <input type="checkbox" checked> thickEnd</label><br/>' +
								'<label> <input type="checkbox" checked> pubMedId</label>' +
							'</div>' +
						'</div>' +
					'</div>' +
					'<div class="col-md-3">' +
						'<div class="panel panel-default">' +
							'<div class="panel-heading">' +
								'<h3 class="panel-title"><label> <input type="checkbox"> DGVa</label></h3>' +
							'</div>' +
							'<div class="panel-body ann-source-box">' +
								'<label> <input type="checkbox"> col1</label> <br/>' +
								'<label> <input type="checkbox"> col2</label><br/>' +
								'<label> <input type="checkbox"> col3</label><br/>' +
								'<label> <input type="checkbox"> col4</label><br/>' +
								'<label> <input type="checkbox"> col5</label><br/>' +
								'<label> <input type="checkbox"> col6</label>' +
							'</div>' +
						'</div>' +
					'</div>' +
					'<div class="col-md-3">' +
						'<div class="panel panel-default">' +
							'<div class="panel-heading">' +
								'<h3 class="panel-title"><label> <input type="checkbox"> Autism_ORoark_Trios</label></h3>' +
							'</div>' +
							'<div class="panel-body ann-source-box">' +
								'<label> <input type="checkbox"> col1</label> <br/>' +
								'<label> <input type="checkbox"> col2</label><br/>' +
								'<label> <input type="checkbox"> col3</label><br/>' +
								'<label> <input type="checkbox"> col4</label><br/>' +
								'<label> <input type="checkbox"> col5</label><br/>' +
								'<label> <input type="checkbox"> col6</label>' +
							'</div>' +
						'</div>' +
					'</div>' +
				'</div>' +
			'</div>' +
		'</div>' +
	'</div>' +
	'<div class="panel panel-default">' +
		'<div class="panel-heading">' +
			'<h3 class="panel-title"><a class="accordion-toggle" data-toggle="collapse" data-parent="#accordion" href="#collapse3">Genome</a></h3>' +
		'</div>' +
		'<div id="collapse3" class="panel-collapse collapse">' +
			'<div class="panel-body">' +
				'<div class="row">' +
					'<div class="col-md-3">' +
						'<div class="panel panel-default">' +
							'<div class="panel-heading">' +
								'<h3 class="panel-title"><label> <input type="checkbox" checked> RefSeq Genes</label></h3>' +
							'</div>' +
							'<div class="panel-body ann-source-box">' +
								'<label> <input type="checkbox" checked> chrom</label> <br/>' +
								'<label> <input type="checkbox" checked> chromStart</label><br/>' +
								'<label> <input type="checkbox" checked> chromEnd</label><br/>' +
								'<label> <input type="checkbox" checked> score</label><br/>' +
								'<label> <input type="checkbox" checked> strand</label><br/>' +
								'<label> <input type="checkbox" checked> txStart</label><br/>' +
								'<label> <input type="checkbox" checked> txEnd</label><br/>' +
								'<label> <input type="checkbox" checked> cdsStart</label><br/>' +
								'<label> <input type="checkbox" checked> cdsEnd</label>' +
							'</div>' +
						'</div>' +
					'</div>' +
					'<div class="col-md-3">' +
						'<div class="panel panel-default">' +
							'<div class="panel-heading">' +
								'<h3 class="panel-title"><label> <input type="checkbox" checked> Ensembl Genes</label></h3>' +
							'</div>' +
							'<div class="panel-body ann-source-box">' +
								'<label> <input type="checkbox" checked> chrom</label> <br/>' +
								'<label> <input type="checkbox" checked> chromStart</label><br/>' +
								'<label> <input type="checkbox" checked> chromEnd</label><br/>' +
								'<label> <input type="checkbox" checked> score</label><br/>' +
								'<label> <input type="checkbox" checked> strand</label><br/>' +
								'<label> <input type="checkbox" checked> txStart</label><br/>' +
								'<label> <input type="checkbox" checked> txEnd</label><br/>' +
								'<label> <input type="checkbox" checked> cdsStart</label><br/>' +
								'<label> <input type="checkbox" checked> cdsEnd</label>' +
							'</div>' +
						'</div>' +
					'</div>' +
					'<div class="col-md-3">' +
						'<div class="panel panel-default">' +
							'<div class="panel-heading">' +
								'<h3 class="panel-title"><label> <input type="checkbox"> GERP</label></h3>' +
							'</div>' +
							'<div class="panel-body ann-source-box">' +
								'<label> <input type="checkbox"> col1</label> <br/>' +
								'<label> <input type="checkbox"> col2</label><br/>' +
								'<label> <input type="checkbox"> col3</label><br/>' +
								'<label> <input type="checkbox"> col4</label><br/>' +
								'<label> <input type="checkbox"> col5</label><br/>' +
								'<label> <input type="checkbox"> col6</label>' +
							'</div>' +
						'</div>' +
					'</div>' +
				'</div>' +
			'</div>' +
		'</div>' +
	'</div>' +
	'<div class="panel panel-default">' +
		'<div class="panel-heading">' +
			'<h3 class="panel-title"><a class="accordion-toggle" data-toggle="collapse" data-parent="#accordion" href="#collapse4">Impact prediction</a></h3>' +
		'</div>' +
		'<div id="collapse4" class="panel-collapse collapse">' +
			'<div class="panel-body">' +
				'<div class="row">' +
					'<div class="col-md-3"><label> <input type="checkbox" checked> SnpEff</label></div>' +
					'<div class="col-md-3"><label> <input type="checkbox"> dbNSFP </label></div>' +
				'</div>' +
			'</div>' +
		'</div>' +
	'</div>' +
	'<div class="panel panel-default">' +
		'<div class="panel-heading">' +
			'<h3 class="panel-title"><a class="accordion-toggle" data-toggle="collapse" data-parent="#accordion" href="#collapse5">Genotype - Phenotype</a></h3>' +
		'</div>' +
		'<div id="collapse5" class="panel-collapse collapse">' +
			'<div class="panel-body">' +
				'<div class="row">' +
					'<div class="col-md-3"><label> <input type="checkbox"> OMIM</label></div>' +
					'<div class="col-md-3"><label> <input type="checkbox"> HGMD </label></div>' +
					'<div class="col-md-3"><label> <input type="checkbox"> ClinVar </label></div>' +
				'</div>' +
			'</div>' +
		'</div>' +
	'</div>' +
	'<div class="panel panel-default">' +
		'<div class="panel-heading">' +
			'<h3 class="panel-title"><a class="accordion-toggle" data-toggle="collapse" data-parent="#accordion" href="#collapse6">Protein networks</a></h3>' +
		'</div>' +
		'<div id="collapse6" class="panel-collapse collapse">' +
			'<div class="panel-body">' +
				'<div class="row">' +
					'<div class="col-md-3"><label> <input type="checkbox"> KEGG</label></div>' +
					'<div class="col-md-3"><label> <input type="checkbox"> HPRD </label></div>' +
					'<div class="col-md-3"><label> <input type="checkbox"> String </label></div>' +
				'</div>' +
			'</div>' +
		'</div>' +
	'</div>' +
	'<div class="panel panel-default">' +
		'<div class="panel-heading">' +
			'<h3 class="panel-title"><a class="accordion-toggle" data-toggle="collapse" data-parent="#accordion" href="#collapse7">Cancer</a></h3>' +
		'</div>' +
		'<div id="collapse7" class="panel-collapse collapse">' +
			'<div class="panel-body">' +
				'<div class="row">' +
					'<div class="col-md-3"><label> <input type="checkbox"> COSMIC</label></div>' +
					'<div class="col-md-3"><label> <input type="checkbox"> TCGA </label></div>' +
					'<div class="col-md-3"><label> <input type="checkbox"> ICGC </label></div>' +
				'</div>' +
			'</div>' +
		'</div>' +
	'</div>' +
	'<div class="panel panel-default">' +
		'<div class="panel-heading">' +
			'<h3 class="panel-title">Pharmacogenomic</h3>' +
		'</div>' +
		'<div id="collapse7" class="panel-collapse collapse">' +
			'<div class="panel-body">' +
				'<div class="row">' +
					'<div class="col-md-3"><label> <input type="checkbox"> DrugBank</label></div>' +
				'</div>' +
			'</div>' +
		'</div>' +
	'</div>' +
'</div>';