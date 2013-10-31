var app = null;

(function () {
    'use strict';

    app = angular.module('SandBox', ['ngResource', 'ui.router', 'ngCookies']);

    app.config(['$routeProvider', function ($routeProvider) {
        var access = routingConfig.accessLevels;

        $routeProvider
            .when('/default', { templateUrl: 'views/index.html', controller: "CommonCtrl", access: access.user })
            .when('/signin', { templateUrl: 'views/signin.html', controller: "SigninCtrl", access: access.anon })
            .when('/signup', { templateUrl: 'views/signup.html', controller: "SignupCtrl", access: access.anon })
            .otherwise({ redirectTo: '/signin' });

    }])
    .run(['$rootScope', '$location', 'authenticationService', function ($rootScope, $location, authenticationService) {
        $rootScope.$on("$routeChangeStart", function (event, next, current) {
            $rootScope.error = null;

            authenticationService.isLogged(function () {
                if (next.access == routingConfig.accessLevels.user) {
                    return;
                }
                $location.path('/default');
            }, function () {
                if (next.access == routingConfig.accessLevels.anon) {
                    return;
                }
                $location.path('/signin');
            });
        });
    }]);
})();