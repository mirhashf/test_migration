(function () {
    'use strict';

    app.controller('CommonCtrl', ['$scope', '$location', 'authenticationService', function ($scope, $location, authenticationService) {
        $scope.data = "test";
        
        $scope.logout = function () {
            authenticationService.logout(function () {
                $location.path("/signin");
            });
        };
    }]);
})();