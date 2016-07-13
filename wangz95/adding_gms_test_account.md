# Adding GMS Test Account

## Steps

### Add an entry for new test account

Add an entry for new test account in
`portal/frontend2/src/test/java/com/bina/seqalto/portal/frontend2/util/UserGenerator.java`:

```java
entries.add("20057, jenkins, [EMAIL], Jenkins CI");
```

### Run `UserGenerator.main()`

New test account will be added to `portal/frontend2/src/main/resources/testData/01-users.json`.

### Commit the changes

Commit the changes introduced by previous two steps, and create a pull request using this commit.

## References
[Pull request 8092](https://github.com/BinaTechnologies/seqalto/pull/8092)

