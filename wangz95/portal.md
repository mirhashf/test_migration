# Portal

## SDK

This module is a thin intelligent layer based on various use cases on top of `portal/api2` REST API provided by the
Portal frontend (PFE).

PortalDataClient
```java
protected final ProxyApi dataResourceApi = new ProxyApi();

public PortalDataClient(URI portalURI, String jwtToken) {
  super(portalURI, jwtToken);
  dataResourceApi.setApiClient(apiClient);
}
```

PortalComputeClient

## Frontend2

## Backend2

### Package `com.bina.seqalto.portal.backend2.asuser`

#### Command.java

#### PortalAsUserService.java

```java
protected static final String DIRECTORY_GLOBAL_PORTALBACKEND = “/usr/lib/bina/portal-backend”;
```

