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

def query(cookies, data_set_id, expression_tree)
  data = JSON.generate({:dataSetId => data_set_id, :expressionTree => expression_tree})
  api_post_json("query?run=true", data, {'Cookie' => cookies})                                      
end

def get_query(id, cookies)
  return api_get("query/#{id}", {'Cookie' => cookies})
end

begin
  cookies = login
  et_proband_gene_tpte = {:subject => "PROBAND",
                           :dataSource => "refgene", 
                           :column => "name2", 
                           :operator => "EQ", 
                           :operand => "TPTE"};
                           
  # run the above query        
  query_response = query(cookies, 1000114, et_proband_gene_tpte)
  
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
rescue Exception => ex
  puts ex.message
end