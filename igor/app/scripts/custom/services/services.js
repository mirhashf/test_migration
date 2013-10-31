(function () {
    'use strict';

    app.factory('messegesService', function () {
        return {
            getMessage: function (type, status) {
                var messages = this[type].filter(function (item) {
                    return item.status == status;
                });

                return (messages[0] && messages[0].message) || "";
            },
            register: [{
                status: 422,
                message: "Invalid email or password does not match criteria"
            }, {
                status: 409,
                message: "email already exists"
            }],
            login: [{
                status: 401,
                message: "user or password is wrong"
            }]
        };
    });

    app.factory('authenticationService', ['$http', '$cookieStore', '$location', function ($http, $cookieStore, $location) {
        return {
            isLogged: function (fnSuccess, fnError) {
                $http({ method: 'GET', url: 'http://ec2-54-212-37-121.us-west-2.compute.amazonaws.com:8080/api/v1/session/user' })
                    .success(function (data, status, headers, config) {
                        fnSuccess && fnSuccess(data, status, headers, config);
                        currentUser = true;
                    })
                    .error(function (data, status, headers, config) {
                        fnError && fnError(data, status, headers, config);
                        currentUser = false;
                    });
            },

            login: function (dataJson, fnError) {
                $http({ method: 'POST', url: 'http://ec2-54-212-37-121.us-west-2.compute.amazonaws.com:8080/api/v1/session', data: dataJson })
                    .success(function (data, status, headers, config) {
                        $location.path("/default");
                    })
                    .error(function (data, status, headers, config) {
                        fnError && fnError(data, status, headers, config);
                    });
            },

            register: function (dataJson, fnError) {
                $http({ method: 'POST', url: 'http://ec2-54-212-37-121.us-west-2.compute.amazonaws.com:8080/api/v1/user', data: dataJson })
                    .success(function (data, status, headers, config) {
                        $location.path("/signin");
                    })
                    .error(function (data, status, headers, config) {
                        fnError && fnError(data, status, headers, config);
                    });
            },

            logout: function (fnSuccess, fnError) {
                $http({ method: 'DELETE', url: 'http://ec2-54-212-37-121.us-west-2.compute.amazonaws.com:8080/api/v1/session', data: '', headers: { "Content-Type": 'application/json' } })
                    .success(function (data, status, headers, config) {
                        $location.path("/signin");
                    })
                    .error(function (data, status, headers, config) {
                        fnError && fnError(data, status, headers, config);
                    });
            }
        };
    }]);
})();