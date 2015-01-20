require 'net/http'
require 'openssl'
require 'json'

BASE_URL = "https://staging.binacloud.com/api/v1"
EMAIL_ADDRESS = "mc@celebs.com" # specific to the api user
PASSWORD = "MileyCyrus1" # specific to the api user
CONTENT_TYPE_JSON = "application/json" 

def get_uri(path)
  URI.parse("#{BASE_URL}/#{path}")
end

def get_https(uri)
  http = Net::HTTP.new(uri.host, uri.port)
  http.use_ssl = true
  http.verify_mode = OpenSSL::SSL::VERIFY_NONE
  
  return http
end

def api_get(path, headers)
  uri = get_uri(path)
  http = get_https(uri)
  return http.get(uri, headers)
end

def api_post(path, data, headers)
  uri = get_uri(path)
  http = get_https(uri)
  return http.post(uri, data, headers)
end

def api_post_json(path, data, headers)
  headers['Content-Type'] = CONTENT_TYPE_JSON
  return api_post(path, data, headers)
end

def login 
  data = JSON.generate({:emailAddress => EMAIL_ADDRESS, :password => PASSWORD})
  response = api_post_json("session", data, {})
  
  raise "Login failed" if response.code != "201"
  
  cookies_in_response = response.get_fields('set-cookie')
  cookies_array = Array.new
  cookies_in_response.each { | cookie |
    cookies_array.push(cookie.split('; ')[0])
  }
  return cookies_array.join('; ')
end

def get_sample_ids_in_project(project_id, cookies)
  response =  api_get("project/#{project_id}", {'Cookie' => cookies})
  project = JSON.parse(response.body);
  groups = project['groups']
  sample_ids = []
  groups.each do |group| 
    samples = group['samples']
    samples.each { |sample| sample_ids.push(sample['id']) }
  end
  sample_ids
end

def dataset_query(cookies, data_set_id, expression_tree)
  data = JSON.generate({:dataSetId => data_set_id, :expressionTree => expression_tree})
  api_post_json("query?run=true", data, {'Cookie' => cookies})                                      
end

def multi_sample_query_with_gene_panel(cookies, sample_ids, expression_tree, refseq_gene_panel)
  group_by = {:dataSource => "refgene", :column => "name2", :values => refseq_gene_panel}
  data = JSON.generate({:sampleIds => sample_ids, :expressionTree => expression_tree, :groupBy => group_by})
  JSON.parse(api_post_json("query?run=true", data, {'Cookie' => cookies}).body)                                      
end

def lollipop_query(cookies, query_id, gene)
  response =  api_get("query/#{query_id}/llpp?gene=#{gene}", {'Cookie' => cookies})
  JSON.parse(response.body)
end

def get_query(id, cookies)
  return api_get("query/#{id}", {'Cookie' => cookies})
end

begin
  cookies = login    
=begin
  et_proband_gene_tpte = {:subject => "PROBAND",
                          :dataSource => "refgene", 
                          :column => "name2", 
                          :operator => "EQ", 
                          :operand => "TPTE"};
                           
  # run the above query    
  
  puts "Running a data set query"
  query_response = dataset_query(cookies, 1000114, et_proband_gene_tpte)
  if query_response.code == '200'
    result = JSON.parse(query_response.body)
    # print total records
    puts "Total records: #{result['totalRecords']}" 
    
    # print the actual variants
    puts "Variants: #{result['variants']}"
    
    # response also contains the query id that was just generated
    query_id = result['queryId']
    puts "Query id: #{query_id}"
    
    # query object can be retrieved by its id
    get_query_by_id_response = get_query(query_id, cookies)
    puts "Query #{query_id}: #{JSON.parse(get_query_by_id_response.body)}"
  else
    raise "Query failed"
  end
=end
  
  # run a multi sample query over all samples in the project
  project_id = 70
  puts "Running a multi sample query #{project_id}"
  sample_ids = get_sample_ids_in_project(project_id, cookies)
  puts "Sample ids in project #{project_id}: #{sample_ids}"
  
  # Run this query for all samples in the project
  multi_sample_exp_tree = {:dataSource => "vcf", :column => "qual", :operator => "GTE", :operand => "50"}
  
  # Group the query results by these genes
  refseq_gene_panel = ["BRCA1", "BRCA2", "BRAF", "KRAS"]
  
  query_response = multi_sample_query_with_gene_panel(cookies, sample_ids, multi_sample_exp_tree, refseq_gene_panel)
  query_id = query_response['queryId']
  
  puts "Multi sample query id: #{query_id}"
  puts "Multi sample num variants map: #{query_response['numVariantsMap']}"
  puts "Multi sample group by results: #{query_response['groupByResults']}"
  
  brca1 = "BRCA1"
  lollipop_response = lollipop_query(cookies, query_id, brca1)
  # lollipop mutation data is per isoform and per amino acid position
  lollipop_mutation_data = lollipop_response['lollipopMutationData']
  puts "Lollipop variants on gene #{brca1} for query #{query_id}: #{lollipop_mutation_data}"
  
  # print the isoform data for BRCA1 in this case
  lollipop_isoform_data = lollipop_response['lollipopIsoformData']
  puts "Lollipop isoforms of gene #{brca1}: #{lollipop_isoform_data}"
  
rescue Exception => ex
  puts ex.message
end