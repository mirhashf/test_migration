(function () {
    'use strict';

    app.controller('SignupCtrl', ['$scope', '$location', 'authenticationService', 'messegesService', function ($scope, $location, authenticationService, messegesService) {
        $scope.redirectToIndex = function () {
            $location.path("/default");
        };
        $scope.redirectToSignin = function () {
            $location.path("/signin");
        };
        $scope.isShow = false;
        $scope.createUser = function () {
            if ($scope.password == $scope.confirmPassword) {
                $scope.isShow = false;
                
                var dataJson = JSON.stringify({
                    emailAddress: $scope.email,
                    password: $scope.password,
                    firstName: "Miley",
                    lastName: "Cyrus",
                    organizationId: 1,
                    occupation: "singer"
                });
                
                authenticationService.register(dataJson, function (data, status, headers, config) {
                    $scope.isShow = true;
                    $scope.message = messegesService.getMessage('register', status);;
                });
            } else {
                $scope.isShow = true;
                $scope.message = "These passwords don't match. Try again?";
            }
        };
    }]);
})();