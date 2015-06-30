#!/usr/bin/python

import sys, os
import boto, boto.ses
import getpass
from optparse import OptionParser
from urlparse import urlparse

def main():

    parser = OptionParser()

    parser.add_option('--zone', default="bina.com", help='Zone', metavar='<domain>', type='string', dest='dns_zone')
    parser.add_option('--server', default="eval.bina.com", help='Eval server', metavar='<url>', type='string', dest='eval_server')
    parser.add_option('--ttl', default=60, help='DNS record TTL', metavar='<seconds>', type='int', dest='dns_ttl')

    (opts, args) = parser.parse_args()

    if (not isinstance(args, (list))) or (not len(args) == 1):
        print ''
        print 'Usage: activate_eval [options] <customer>'
        print ''
        sys.exit()

    customer = args[0]

    print "using customer: {0}".format(customer)
    print "using TTL: {0}".format(opts.dns_ttl)
    print "using zone: {0}".format(opts.dns_zone)

    if (raw_input('Proceed using these parameters? [Y/n]: ') != 'Y'):
        sys.exit()

    r53_region = 'us-east-1'

    r53_connection = boto.route53.connect_to_region(r53_region)
    dns_zone = r53_connection.get_zone(opts.dns_zone)

    dns_zone.add_record("CNAME", "{0}.{1}".format(customer, opts.dns_zone), opts.eval_server)


if __name__ == '__main__':
    main()

