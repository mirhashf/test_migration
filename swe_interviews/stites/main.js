var main = angular.module('main', []);

main.controller('resCtrl', function($scope, $http, $q){
	$http({ method: 'GET',
		    url:'http://localhost:9000/api/mutations?pagesize=32'})
	.then(function(res){
		$scope.response = res.data.objects;
		for(var i=0; i < $scope.response.length; i++){
			(function(i){
				$http({method: 'GET',
					   url: 'http://localhost:9000/api/mutations/'+$scope.response[i].id+'/annotations'})
				.then(function(res){
					$scope.response[i].annotations = res.data;
				})
			})(i)
		}
	})
});