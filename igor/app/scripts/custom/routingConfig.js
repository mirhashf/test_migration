(function (exports) {

    var userRoles = {
        publicRole: 1, // 001
        userRole: 2, // 010
        adminRole: 4 // 100
    };

    exports.userRoles = userRoles;
    exports.accessLevels = {
        publicLevel: userRoles.publicRole | // 111
            userRoles.userRole |
            userRoles.adminRole,
        anonLevel: userRoles.publicRole, // 001
        userLevel: userRoles.userRole | // 110
            userRoles.adminRole,
        adminLevel: userRoles.adminRole // 100
    };
})(typeof exports === 'undefined' ? this['routingConfig'] = {} : exports);