function getNewFilterRow(rowId, xId, isFirst) {
	return '<div class="row" style="margin-bottom: 3px" id="'
			+ rowId
			+ '">'
			+ '<div class="col-md-12">'
			+ '<label'
			+ (isFirst ? "" : ' style="visibility: hidden"')
			+ '>Filter</label> '
			+ '<select>'
			+ '<option>---</option>'
			+ '<optgroup label="Small Variations">'
			+ '<option>dbSNP</option>'
			+ '<option>1000 Genomes</option>'
			+ '<option>HapMap</option>'
			+ '</optgroup>'
			+ '<optgroup label="Structural Variations">'
			+ '<option>dbVar</option>'
			+ '<option>DGV</option>'
			+ '<option>DGVa</option>'
			+ '<option>Autism_ORoark_Trios</option>'
			+ '</optgroup>'
			+ '<optgroup label="Genome">'
			+ '<option>Ensembl</option>'
			+ '<option>RefSeq</option>'
			+ '<option>GERP</option>'
			+ '</optgroup>'
			+ '<optgroup label="Impact Prediction">'
			+ '<option>SnpEff</option>'
			+ '<option>dbNSFP</option>'
			+ '</optgroup>'
			+ '<optgroup label="Genotype - Phenotype">'
			+ '<option>OMIM</option>'
			+ '<option>HGMD</option>'
			+ '<option>ClinVar</option>'
			+ '</optgroup>'
			+ '<optgroup label="Protein Networks">'
			+ '<option>KEGG</option>'
			+ '<option>HPRD</option>'
			+ '<option>String</option>'
			+ '</optgroup>'
			+ '<optgroup label="Cancer">'
			+ '<option>COSMIC</option>'
			+ '<option>TCGA</option>'
			+ '<option>ICGC</option>'
			+ '</optgroup>'
			+ '<optgroup label="Pharmacogenomic">'
			+ '<option>DrugBank</option>'
			+ '</optgroup>'
			+ '</select> '
			+ '<select>'
			+ '<option>---</option>'
			+ '<option>name</option>'
			+ '<option>chrom</option>'
			+ '<option>chromStart</option>'
			+ '<option>chromEnd</option>'
			+ '<option>class</option>'
			+ '<option>func</option>'
			+ '</select> '
			+ '<select>'
			+ '<option>---</option>'
			+ '<option>=</option>'
			+ '<option>!=</option>'
			+ '<option>&gt;</option>'
			+ '<option>&gt;=</option>'
			+ '<option>&lt;</option>'
			+ '<option>&lt;=</option>'
			+ '<option>IN</option>'
			+ '<option>NOT IN</option>'
			+ '<option>LIKE</option>'
			+ '</select> '
			+ '<input type="text" size="80"> '
			+ '<a href="" id="'
			+ xId
			+ '"><span class="glyphicon glyphicon-remove"></span></a> '
			+ (isFirst ? '<button class="btn btn-primary" id="add-new-filter">Add</button>'
					: "") + '</div>' + '</div>';
}

var newFilterRowCounter = 0;

function addNewRow(panel, isFirst) {
	newFilterRowCounter++;
	var rowId = 'filter-row-' + newFilterRowCounter;
	var xId = 'remove-row-' + newFilterRowCounter;
	panel.append(getNewFilterRow(rowId, xId, isFirst));
	var rowElem = $('#' + rowId);
	var xElem = $('#' + xId);
	xElem.click(function() {
		rowElem.remove();
		return false;
	});
}