# Release Tests

## Installing Lombok plugin for Eclipse

[Instructions](http://jnb.ociweb.com/jnb/jnbJan2010.html)

## Running Release Tests
In order to run the release tests, you'll generally need to update several local boxes, launch and tests and then
monitor them.
While this can all be done with a single command, there are many pieces of software involved.
All of these steps should be run from **t-rex**, and in a **tmux** session as they will run for many hours and should
not be interrupted.

## Steps

### Reserve a Docker compute resource

Reserve a Docker compute resource (say `tehran-02-06`) on [Compute Resources](https://binatechnologies.atlassian.net/wiki/display/ENG/Compute+Resources).

### Set up SSH tunneling from laptop to **t-rex**
```bash
tmux

ssh -L 8081:localhost:8081 t-rex -A
```

### Tear down previous active ITL deployment

Check if there is an active ITL deployment running:
```bash
ssh t-rex -A

ssh build-00

ps -ef|grep wangz95
```

Tear down it from your laptop:
```bash
java -cp test/integration/target/assembly/test-integration-launcher.jar com.bina.seqalto.test.integration.IntegrationTest_Teardown --rigtype "BinaDevelopment" --rigargs "--gitURL git@github.com:BinaTechnologies/seqalto.git develop tehran-01-06" --username wangz95 --noprompt
```

### Deploy ITL to the compute resource you reserved

Deploy ITL from your laptop
```bash
tmux

java -jar test/integration/target/assembly/test-integration-launcher.jar --rigtype "BinaDevelopment" --rigargs "--gitURL git@github.com:BinaTechnologies/seqalto.git develop tehran-02-06" --username wangz95 --noprompt --keepup
```

Write down the URL from following line:
```
2016-05-23 17:14:17,719 INFO  c.b.s.t.i.tester.IntegrationTester - Setup client tunnel as http://localhost:18973
```

**NOTE:** If you disconnected from above session, you need to set up port forwarding using `ssh -D 7777 t-rex`, and use
`http://build-00:18973` to access the portal site in browser from your laptop.

### Get JWT token

#### Using Web browser

Install the "Postman - REST Client shortcut" extension in your Chrome.

Log in by sending a POST request with your credentials to above URL plus `/api/session`, i.e.
`http://localhost:18973/api/session` (or `http://build-00:18973/api/session`) in this case:
```
{"name":"jerryw","password":"b"}
```

Send a POST request with your credentials to above URL plus `/api/jwt`, i.e. `http://localhost:18973/api/jwt`
(or `http://build-00:18973/api/jwt`) in this case:
```
{"name":"jerryw","password":"b"}
```

Write down the JWT token returned in response.

#### Using curl

You can do that using `curl` and `python` as following:
```bash
curl -s -H "Content-Type: application/json" -X POST -d '{"name":"jerryw","password":"b"}' http://localhost:18973/api/session -c cookieFile

curl -s -H "Content-Type: application/json" -X POST -d '{"name":"jerryw","password":"b"}' http://localhost:18973/api/jwt -b cookieFile | python -c 'import json,sys;resp=json.load(sys.stdin);print resp["data"]["token"]'
```

### Running test case in Eclipse

Set up a SSH tunnel like following:
```bash
ssh -tt -L 18973:localhost:18973 t-rex -A ssh -L 18973:localhost:18973 build-00
```
**NOTE:** Port number needs to be updated in above command according to your frontend port number.

In Eclipse, run `TCx40_0.java` with **Run As > Run Configurations...** configured with following system properties in
**VM arguments** in **Arguments** tab:
```
-DENV_URL=http://localhost:18973 -DJWT_TOKEN=eyJhbGciOiJIUzUxMiJ9.eyJzdWIiOiJqZXJyeXciLCJleHAiOjE0NjY2MzA4NzcsImlzcyI6IkJpbmEgR01TLVJBVkUiLCJpYXQiOjE0NjQwMzg4NzcsImp0aSI6IkFQRFpaSnhSU254Qms0QnBaa3Z6UHJISiIsIkZPUkNFX1JPTEUiOmZhbHNlfQ.e5bCqeUz5uqOHHmQ9Cy0vIh2cG6cZECUiVDyFVD5zmvhwyrtR99oh6XrnTuvZ2CV1aRZ8HzcocHaRxK-jdQJbw -DTEST_USER_EMAIL=[YOUR EMAIL] -DTEST_LOCAL_DIR=releasetest
```

As you can see, following system properties are required in above **Run Configuration**:
* ENV_URL
* JWT_TOKEN
* TEST_USER_EMAIL
* TEST_LOCAL_DIR

### Running test case in IntelliJ IDEA

Set up a SSH tunnel like following:
```bash
ssh -tt -L 19013:localhost:19013 t-rex -A ssh -L 19013:localhost:19013 lima-00
```
**NOTE:** The port number can be figured out by examining the `portal-frontend` process running on `lima-00`.

In IntelliJ IDEA, right click `TCx51.java`, click `Create 'TCx51'...`.
Add following system properties to **VM options** in **Configuration** tab:
```
-DENV_URL=https://localhost:19013 -DJWT_TOKEN=eyJhbGciOiJIUzUxMiJ9.eyJzdWIiOiJqZXJyeXciLCJleHAiOjE0Njk2NDA2NTQsImlzcyI6IkJpbmEgR01TLVJBVkUiLCJpYXQiOjE0NjcwNDg2NTQsImp0aSI6IkFNTTBFUE83UjNoSGlZTHpnRGpWb21GWCIsIkZPUkNFX1JPTEUiOmZhbHNlfQ.PvO7krAVzzi626b1c9zzT44T_mvk5obtafqnJB-9fMWCHQkDWqe0J042O_cX9A0VVd3rGbCtN_sMyCYHgluHoA -DTEST_USER_EMAIL=[YOUR EMAIL] -DTEST_LOCAL_DIR=releasetest
```
Click "OK" button to save the configuration. Now you can right click `TCx51.java` then click `Run 'TCx51'` to run it.

### Running test suite from command line

Set up a SSH tunnel like following:
```bash
ssh -tt -L 18973:localhost:18973 t-rex -A ssh -L 18973:localhost:18973 build-00
```

```bash
java -cp ./test-release-3.3.0-SNAPSHOT.jar:lib/loomis2libs/* -DENV_URL=https://lima-00 org.junit.runner.JUnitCore com.bina.seqalto.test.release.dna.TCx40_8
```

Run `mvn test surefire-report:report` with following system properties:
```bash
mvn -DENV_URL=http://localhost:18973 -DJWT_TOKEN=[JWT TOKEN] -DTEST_USER_EMAIL=[YOUR EMAIL] -DTEST_LOCAL_DIR=/tmp/bina_test_suite -Dtest=AllFastTestsSuite test surefire-report:report
```

### Polling for state update after job submission

```bash
mvn -DENV_URL=http://localhost:18973 -DJWT_TOKEN=[JWT TOKEN] -DTEST_USER_EMAIL=[YOUR EMAIL] -DTEST_LOCAL_DIR=/tmp/bina_test_suite -Dtest=TCTestStateUpdater test surefire-report:report
```

### Job validation

After `TCTestStateUpdater` completed, run `AllFastTestsValidationSuite` to do the job validation.
```bash
mvn -DENV_URL=http://localhost:18973 -DJWT_TOKEN=[JWT TOKEN] -DTEST_USER_EMAIL=[YOUR EMAIL] -DTEST_LOCAL_DIR=/tmp/bina_test_suite -Dtest=AllFastTestsValidationSuite test surefire-report:report
```


## Test case naming convention
* `TCx40` for DNA
* `TCx50` for RNA
* `TCx80` for Somatic
* `TCx90` for Multisample

## Old Release Tests

All old release tests reside in
[seqalto/test/integration/src/test/json](https://github.com/BinaTechnologies/seqalto/tree/develop/test/integration/src/test/json).

## Test plan

Take a look at `test/integration/src/test/json/workgen-dna/TCx40.0.js`.

Result monitoring

## Code Reading

```java
io.swagger.client.model.PortalJob.StatusEnum.SUCCESSFUL;

io.swagger.client.model.ASecondClassApiobject
```
Generated code resides in `portal/api2/src/main/java/io/swagger/client/model`

```java
public enum StatusEnum {
  UNKNOWN("Unknown"),
  QUEUED("Queued"),
  RUNNING("Running"),
  WAITING("Waiting"),
  KILLED("Killed"),
  SUCCESSFUL("Successful"),
  ERROR("Error");

  private String value;

  StatusEnum(String value) { this.value = value; }
}
```

```java
com.fasterxml.jackson.databind.JsonNode
com.fasterxml.jackson.databind.ObjectMapper
com.fasterxml.jackson.databind.SerializationFeature
com.google.common.base.Supplier
```

### Validation

#### TestState

Defined in `test/release/src/test/java/com/bina/seqalto/test/release/validation/TestState.java`
```java
public class TestState {
  public TestState() {
  }
}
```

#### ValidationState

Defined in `test/release/src/test/java/com/bina/seqalto/test/release/validation/ValidationState.java`
```java
public class ValidationState {
  protected final LinkedHashSet<TestState> testStates = new LinkedHashSet<>();

  /**
   * Map between Class:TestMethodName and TestState.
   */
  protected final HashMap<Tuple2<Class, String>, TestState> testMethodState = new HashMap<>();
}
```

### Test

#### TCTest

The directory where all the configuration JSON files reside in
```java
public static final String CONFIGURATION_BASE_PATH = "/com/bina/seqalto/test/release/configuration/";
```

`test/release/src/test/java/com/bina/seqalto/test/release/rna/RNAFastTestSuite.java`
```java
@RunWith(Suite.class)
@SuiteClasses({TCx50.class, TCx51.class})
public class RNAFastTestSuite {
}
```

```java
test50dot0

createNFSOutputWorkflow(System system) {
  final Workflow workflow = new Workflow();
  workflow.setDebug(new Debug(Debug.SCRATCH.TRUE));
  workflow.setMounts(createMounts());
  workflow.setSystem(system);
  setUserAndNotification(workflow);
  createWorkflowOutput(workflow, "755", "mount:river/output", "Parent-Project-000000000");
  workflow.setModules(new Object());
  return workflow;
}
```

```java
/**
* Method to deserialize JSON content as tree expressed using set of JsonNode instances.
*/
public JsonNode readTree(InputStream in) {

}
```

```java
package com.bina.seqalto.qe.workflow.model;

public class System {
  final int workers;
  Machine machine;
}

public class Workflow {
  protected System system;
}
```

Defined `test/release/src/test/java/com/bina/seqalto/test/release/base/TCTest.java`
```java
@Slf4j
public abstract class TCTest {
  // Global validation state
  protected final static ValidationState validationState;

  protected static String getSystemPropertyWithNullCheck(String propertyName) {
  }

  protected static void createWorkflowOutput() {
    final Output output = new Output();
    output.setPermissions(permissions);
    output.setUri(uri);
    output.setProject(portalDataClient.listProject(projectId));
    workflow.setOutput(output);

    if (log.isDebugEnabled()) {
      log.debug(objectMapper.writeValueAsString(output);
    }
  }

  protected static AnalysisApi runAnalysis() {
    final AnalysisApi analysisWithWorkflowCreated = portalComputeClient.submitAnalysis(analysis, analysisCreated.getId());
  }
}
```

```java
protected Boolean expectedStatus(String testMethodName, Class className) {
  final Tuple2<Class, String> testMethodKey = new Tuple2<>(className, testMethodName);
  if (validationState.getTesetMethodState().containsKey(testMethodKey)) {
    TestState testState = validationState.getTestMethodState().get(testMethodKey);
    return testState.getExpectedStatus().equals(testState.getStatus()) ? Boolean.TRUE : Boolean.FALSE;
  }
  return null;
}
```

```java
if (checkExpectedStateIfValidating(() -> name.getMethodName())) return;
```

```java
public boolean checkExpectedStateIfValidating(final Supplier<String> testMethodNameSupplier)
```

```java
TCx93.testAnalysis
TCTest.createAndRunAnalysis
TCTest.runAnalysis
PortalComputeClient.submitAnalysis
APortalSDKClient.invokeWithRetries
PortalComputeCilent$$Lambda$6
PortalComputeClient.lambda$submitAnalysis
io.swagger.client.api.AnalysisactionresouceApi.runUsingPUT (portal/api2/src/main/java/io/swagger/client/api)
io.swagger.client.ApiClient.invokeAPI (portal/api2/src/main/java/io/swagger/client/ApiClient.java)
```

