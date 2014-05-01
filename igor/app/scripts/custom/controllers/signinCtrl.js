(function () {
    'use strict';

    app.controller('SigninCtrl', ['$scope', '$location', 'authenticationService', 'messegesService', function ($scope, $location, authenticationService, messegesService) {
        $scope.redirectToIndex = function () {
            $location.path("/default");
        };
        $scope.redirectToSignup = function () {
            $location.path("/signup");
        };
        $scope.isShow = false;
        $scope.signin = function () {
            var dataJson = JSON.stringify({
                emailAddress: $scope.email,
                password: $scope.password
            });

            authenticationService.login(dataJson, function (data, status, headers, config) {
                $scope.isShow = true;
                $scope.message = messegesService.getMessage('login', status);;
            });
        };
    }]);
})();