#!/usr/bin/env python

import urllib2
import json
from optparse import OptionParser
import sys

class Client(object):
    """Simple REST Client to communicate with Bina Portal.


    Attributes:
        url: A String REST API URL handle
        usr: A String for Bina Portal Username
        passwd: A String for Bina Portal Password
        verbose: A Boolean indicating whether to print detailed logs
    """
    def __init__(self, url, user, passwd, verbose = False):
        self.url = url
        self.user = user
        self.passwd = passwd
        self.token = None
        self.verbose = verbose
        
    def login(self):
        session = self.post({"user": {"name" :self.user, "password": self.passwd}}, "session")
        if session is not None and "id" in session:
            self.token = session["id"]
        elif not "id" in session:
            print "No valid token returned, exiting..."
            exit(1);
    
    def logout(self):
        self.delete ("session/" + self.token)

    def get(self, *args, **kwargs):
        result = self.query('GET', None, *args, **kwargs)
        return result

    def post(self, data, *args, **kwargs):
        result = self.query('POST', data, *args, **kwargs)
        return result

    def put(self, data, *args, **kwargs):
        result = self.query('PUT', data, *args, **kwargs)
        return result

    def delete(self, *args, **kwargs):
        result = self.query('DELETE', None, *args, **kwargs)
        return result

    def query(self, method, data, *args, **kwargs):
        url = '%s/api/%s/' % (self.url, '/'.join(args))
        if len(kwargs):
            url += '?' + '&'.join(['%s=%s' % (k, v) for k, v in kwargs.iteritems()])
        return self._execute(method, url, args, data, kwargs)

    def _execute(self, method, url, args, data, kwargs):
        if method in ['GET', 'DELETE']:
            r = urllib2.Request(url)
        else:
            r = urllib2.Request(url, json.dumps(data))

        r.add_header('Content-Type', 'application/json')
        if(self.token is not None):
            r.add_header('Cookie', "token=" + self.token)
        r.get_method = lambda: method
        o = urllib2.OpenerDirector()
        o.add_handler(urllib2.HTTPHandler())
        o.add_handler(urllib2.HTTPSHandler()) 

        result = o.open(r)
        retdata = result.read()
        retinfo = result.info()
        
        if self.verbose:
            print method, url
            print "DATA", data
            print 'RESULT HEADERS'
            print retinfo
            print 'RESULT DATA'
            print repr(retdata)
            print '%s %s' % (method, url)
        if retdata != "":
            try:
                data = json.loads(retdata)
            except: 
                print retdata
        return data

def check_or_error(var, message, parser):
    if not var:
        print >> sys.stderr, message 
        parser.print_help(sys.stderr)
        exit(1)

def try_pretty_print_json(json_object):
    try:
        return json.dumps(json_object, sort_keys=True, indent=4, separators=(',', ': '))
    except:
        return json_object 
    

def main():
    parser = OptionParser()
    parser.add_option("--get", action="store_const", help="Get request", dest="type", const="get", default="get")
    parser.add_option("--post", action="store_const", help="Post request", dest="type", const="post")
    parser.add_option("--delete", action="store_const", help="Delete request", dest="type", const="delete")
    parser.add_option("--put", action="store_const", help="Put request", dest="type", const="put")
    parser.add_option("--baseurl", type="string", dest="baseurl", help="Base REST URL")
    parser.add_option("--resource", type="string", dest="resource", help="REST resource")
    parser.add_option("--username", type="string", dest="username", help="Username")
    parser.add_option("--password", type="string", dest="password", help="Password")
    parser.add_option("--data", type="string", dest="data_file", help="Data file to be passed to POST or PUT request", metavar="FILE")
    parser.add_option("-v", "--verbose", action="store_true", help="Verbose mode.", dest="verbose", default=False)
    
    (opts, args) = parser.parse_args()
    
    check_or_error(opts.username, "No username specified", parser)
    check_or_error(opts.password, "No username specified", parser)
    check_or_error(opts.baseurl, 'No base REST URL specified.', parser)
    check_or_error(opts.resource, 'No REST resource specified.', parser)
     

    c = Client(opts.baseurl, opts.username, opts.password, opts.verbose)
    c.login()

    def url_query(handler):
        print try_pretty_print_json(handler(opts.resource))

    def data_query(handler):
        data = "{}"
        if opts.data_file:
            with open (opts.data_file, "r") as datafile:
                data = datafile.read().replace('\n', '')
        else:
            data = sys.stdin.read().replace('\n', '')
        print try_pretty_print_json(handler(json.loads(data), opts.resource))
    
    if opts.type == "get":
        url_query(c.get)
        
    if opts.type == "delete":
        url_query(c.delete)
                
    if opts.type == "post":
        data_query(c.post)
        
    if opts.type == "put":
        data_query(c.put)
        
    c.logout()   

   
if __name__ == '__main__':
    main()
