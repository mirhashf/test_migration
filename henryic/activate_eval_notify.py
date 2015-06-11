#!/usr/bin/python

import sys, os
#import re
#import time
#from datetime import datetime
import boto
import getpass
from optparse import OptionParser
from urlparse import urlparse

def main():

    parser = OptionParser()

    parser.add_option('--portal', default='https://eval.bina.com', help='Portal URL', metavar='<URL>', type='string', dest='portal_url')
    parser.add_option('--user', default='bina-admin', help='Administrative portal user', metavar='<username>', type='string', dest='portal_user')
    parser.add_option('--password', help='Administrative portal user password', metavar='<password>', type='string', dest='portal_pass')

    (opts, args) = parser.parse_args()

    if (not isinstance(args, (list))) or (not len(args) == 2):
        print ''
        print 'Usage: activate_eval [options] <customer> <username> <e-mail>'
        print ''
        sys.exit()

    customer = args[0]
    username = args[1]
    cusemail = args[2]

    print opts.portal_url
    print opts.portal_user
    print opts.portal_pass
    print customer
    print username
    print cusemail

    ses_region = 'us-east-1'
    ses_sender = 'mail-noreply@bina.com'

    message_subject = 'Account activated for the Bina RAVE Platform evaluation'
    message_format = 'html'
    message_text_body = 'Hi {0},\n\nYour account for the Bina RAVE evaluation platform has been activated.\n\nYou may now log in at https://{0}.bina.com/#.\n\nPlease do not hesitate to contact us at support@bina.com regarding any concerns.\n\nThe Bina Team'.format(username,customer)
    message_html_body = '<html><body><p>Hi <b>{0}</b>,<br><br>Your account for the Bina RAVE evaluation platform has been activated.<br></p><p>You may now login by clicking this <a href="https://{1}.bina.com/#">link</a>.</p><p>Please do not hesitate to <a href="mailto:support@bina.com">contact us</a> regarding any concerns.</p><p>The Bina Team</p></html>'.format(username,customer)

    ses_connection = boto.ses.connect_to_region(ses_region)

    ses_connection.send_email(source=ses_sender, subject=message_subject, body=None, to_addresses=[cusemail], bcc_addresses=['henryic@bina.com'], format=message_format, text_body=message_text_body, html_body=message_html_body)


if __name__ == '__main__':
    main()

