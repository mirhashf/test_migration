#!/usr/bin/env python

import sys, os
import netconfig
from NetworkInterface import NetworkInterface
from bottle import Bottle, run, get, post, request, route, static_file

@route('/useless_form')
def useless_form():
    current_config = netconfig.get_current_ifaces()
    iface_eth0 = current_config[1]
    return '''<p>
            Current address: %s<br>
            Current netmask: %s<br>
            Current gateway: %s</p>
            <form method="POST" action="/useless_form">
                <p>New address: <input name="address" type="text" /><br>
                   New netmask: <input name="netmask" type="text" /><br>
                   New gateway: <input name="gateway" type="text" /><br>
                   <input type="submit" /></p>
            </form>''' % (iface_eth0.get_address(), iface_eth0.get_netmask(), iface_eth0.get_gateway())

@post('/useless_form')
def useless_submit():
    new_address = request.forms.get('address')
    new_netmask = request.forms.get('netmask')
    new_gateway = request.forms.get('gateway')

    config_command = './netconfig.py --iface=eth0 --address=%s --netmask=%s --gateway=%s' % (new_address, new_netmask, new_gateway)

    os.system(config_command)

    return '''Update submitted. Click <a href="useless_form">here</a> to return.'''

@route('<filepath:path>')
def server_static(filepath):
    return static_file(filepath, root='./static')

run(host='0.0.0.0', port=8080, debug=True)

