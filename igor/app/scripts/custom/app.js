var app = null;

(function () {
    'use strict';

    app = angular.module('SandBox', ['ui.router']).
        config(function ($stateProvider, $urlRouterProvider) {

            // For any unmatched url, send to /route1
            $urlRouterProvider.otherwise("/route1");

            $stateProvider
                .state('index', {
                    templateUrl: 'index.html',
                    controller: function ($scope) {
                        $scope.contacts = [{ name: 'Alice' }, { name: 'Bob' }];
                    }
                })
                .state('index.default', {
                    views: {
                        'list': {
                            templateUrl: 'contacts/list.html',
                            controller: function ($scope) {
                                $scope.contacts = [{ name: 'Alice' }, { name: 'Bob' }];
                            }
                        },
                        'list2': {
                            templateUrl: 'contacts/list2.html'
                        }
                    }
                });

            $stateProvider
                .state('signin', {
                    templateUrl: 'index.html',
                    controller: function ($scope) {
                        $scope.contacts = [{ name: 'Alice' }, { name: 'Bob' }];
                    }
                })
                .state('signin.default', {
                    views: {
                        'list': {
                            templateUrl: 'contacts/list.html',
                            controller: function($scope) {
                                $scope.contacts = [{ name: 'Alice' }, { name: 'Bob' }];
                            }
                        },
                        'list2': {
                            templateUrl: 'contacts/list2.html'
                        }
                    }
                });

            //.state('contacts', {
            //    abstract: true,
            //    templateUrl: 'contacts.html',
            //    controller: function ($scope) {
            //        $scope.contacts = [{ name: 'Alice' }, { name: 'Bob' }];
            //    }
            //})
            //.state('contacts.list', {
            //    url: 'list',
            //    templateUrl: 'contacts/list.html'
            //})
            //.state('contacts.list.details', {
            //    url: 'details',
            //    templateUrl: 'contacts/list.html'
            //});
        });

    app.controller('MainCtrl', function ($state) {
        $state.transitionTo('contacts.list');
    });


    //.config(['$routeProvider', '$locationProvider', '$httpProvider', function ($routeProvider, $locationProvider, $httpProvider) {
    //    var access = routingConfig.accessLevels;

    //    $routeProvider
    //            .when('/signin', { template: 'views/signin.html', controller: "SigninCtrl", access: access.public })
    //            .when('/signup', { template: 'views/signup.html', controller: "SigninCtrl", access: access.anon })
    //            .otherwise({ redirectTo: '/signin' });

    //    //$locationProvider.html5Mode(true);

    //    var interceptor = ['$location', '$q', function ($location, $q) {
    //        function success(response) {
    //            return response;
    //        }

    //        function error(response) {
    //            if (response.status === 401) {
    //                $location.path('/signin');
    //                return $q.reject(response);
    //            }
    //            else {
    //                return $q.reject(response);
    //            }
    //        }

    //        return function (promise) {
    //            return promise.then(success, error);
    //        };
    //    }];

    //    $httpProvider.responseInterceptors.push(interceptor);
    //}])
    //.run(['$rootScope', '$location', 'Auth', function ($rootScope, $location, Auth) {
    //    $rootScope.$on("$routeChangeStart", function (event, next, current) {
    //        $rootScope.error = null;
    //        if (!Auth.authorize(next.access)) {
    //            if (Auth.isLoggedIn()) {
    //                $location.path('/');
    //            } else {
    //                $location.path('/signin');
    //            }
    //        }
    //    });

    //}]);
})();