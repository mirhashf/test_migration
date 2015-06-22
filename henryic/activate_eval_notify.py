#!/usr/bin/python

import sys, os
#import re
#import time
#from datetime import datetime
import boto, boto.ses
import getpass
from optparse import OptionParser
from urlparse import urlparse

def main():

    parser = OptionParser()

    parser.add_option('--portal', default='https://eval.bina.com', help='Portal URL', metavar='<URL>', type='string', dest='portal_url')
    parser.add_option('--user', default='bina-admin', help='Administrative portal user', metavar='<username>', type='string', dest='portal_user')
    parser.add_option('--password', help='Administrative portal user password', metavar='<password>', type='string', dest='portal_pass')

    (opts, args) = parser.parse_args()

    if (not isinstance(args, (list))) or (not len(args) == 3):
        print ''
        print 'Usage: activate_eval [options] <customer> <username> <e-mail>'
        print ''
        sys.exit()

    customer = args[0]
    username = args[1]
    cusemail = args[2]

    print "Using portal URL: {0}".format(opts.portal_url)
    print "Using portal user: {0}".format(opts.portal_user)
    print "using customer: {0}".format(customer)
    print "using username: {0}".format(username)
    print "using customer email: {0}".format(cusemail)

    if (raw_input('Proceed using these parameters? [Y/n]: ') != 'Y'):
        sys.exit()

    ses_region = 'us-east-1'
    ses_sender = 'mail-noreply@bina.com'

    message_subject = 'Account activated for the Bina RAVE Platform evaluation'
    message_format = 'html'
    message_text_body = 'Hi {0},\n\nYour account for the Bina RAVE evaluation platform has been activated.\n\nYou may now log in at https://{0}.bina.com/#. If you were previously logged in, you may need to log out and log back in for the change to take effect.\n\nPlease do not hesitate to contact us at support@bina.com regarding any concerns.\n\nThe Bina Team'.format(username,customer)
    message_html_body = '<html><body><p>Hi <b>{0}</b>,<br><br>Your account for the Bina RAVE evaluation platform has been activated.<br></p><p>You may now login by clicking this <a href="https://{1}.bina.com/#">link</a>. If you were previously logged in, you may need to log out and log back in for the change to take effect.</p><p>Please do not hesitate to <a href="mailto:support@bina.com">contact us</a> regarding any concerns.</p><p>The Bina Team</p></html>'.format(username,customer)

    ses_connection = boto.ses.connect_to_region(ses_region)

    ses_connection.send_email(source=ses_sender, subject=message_subject, body=None, to_addresses=[cusemail], bcc_addresses=['henryic@bina.com'], format=message_format, text_body=message_text_body, html_body=message_html_body)


if __name__ == '__main__':
    main()

