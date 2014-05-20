var app = angular.module('main', ['ngResource', 'ngRoute']);

var apibaseUrl = "http://staging.binacloud.com:9000/api";

app.controller('MainCtrl', function ($scope, $location, $http, $interval) {
    var timer;
    var qr = /^(Query )(\d+)( has been created)$/g
    
    $scope.getNotificationsPeriodically = function () {
        // call this once before the timer is kicked off 
        // otherwise there is an initial wait until the first call
        getNotifications();

        if (timer) {
            $interval.cancel(timer);
        }

        // Refresh every 10 minutes
        timer = $interval(getNotifications, 10*60*1000);
    };

    var getNotifications = function () {
            $scope.loading = true;
            $scope.errorMessage = '';
            $scope.notifications = [];

            $http.get([apibaseUrl, 'activity'].join('/'), {
                params: {
                  targetEmails: $scope.targetEmails
                }
            }).success(function (response, status) {
                $scope.notifications = response;
                $scope.notifications.forEach(function (n) {
                	if (n.message.match(qr)) {
                		n.isQuery = true;
                	}
                })
                $scope.loading = false;
            }).error(function (response, status) {
                $scope.errorMessage = response.errors ? response.errors[0] : 'Error';
                $scope.loading = false;
            });
        };
        
    $scope.getQuery = function (message) {
    	$scope.loading = true;
        $scope.errorMessage = '';
        $scope.query = '';
        
        var queryId = '';
        var groups =  /^(Query )(\d+)( has been created)$/g.exec(message);
        if (groups && groups.length >= 3) {
            queryId = groups[2];
            
            $http.get([apibaseUrl, 'queries', queryId].join('/')).success(function (response, status) {
                $scope.query = response;
                $scope.loading = false;
            }).error(function (response, status) {
                $scope.errorMessage = response.errors ? response.errors[0] : 'Error';
                $scope.loading = false;
            });
        } else {
        	$scope.loading = false;
        }
    };

});

app.filter('formatDate', function () {
    return function (dateStr) {
        return dateStr != null ? moment(dateStr).fromNow() : "";
    };
});