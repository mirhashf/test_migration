// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;

import javax.annotation.PostConstruct;
import javax.annotation.PreDestroy;
import javax.inject.Inject;

import org.apache.zookeeper.CreateMode;
import org.apache.zookeeper.KeeperException;
import org.apache.zookeeper.KeeperException.Code;
import org.apache.zookeeper.WatchedEvent;
import org.apache.zookeeper.Watcher;
import org.apache.zookeeper.Watcher.Event.EventType;
import org.apache.zookeeper.ZooDefs.Ids;
import org.apache.zookeeper.ZooKeeper;
import org.apache.zookeeper.data.ACL;
import org.apache.zookeeper.data.Stat;
import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.Version;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;
import org.codehaus.jackson.map.module.SimpleModule;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Service;

/**
 * A singleton service for each app-server that is the main point of interaction with ZooKeeper.
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 */
@Service
public class ZooKeeperService implements Watcher {
  private static final String SEP = "/";
  private static final String ELECTIONS_NODE_STRING = "elections";
  private static final String RUNNERS_NODE_STRING = "runners";
  private static final String DONE_RUNNERS_NODE_STRING = "done_runners";
  private static final String JOBS_NODE_STRING = "jobs";
  private static final String JOB_NODE_STRING = "job.";
  private static final String NODE_PREFIX_STRING = "node.";

  private static class Config {
    private static boolean isLeader;
    private static String hostname;
    private static ZooKeeper zk;
    private static String boxNode;
    private static String jobsNodePath;
    private static String leaderNodePath;
    private static String myPath;
    private static String electionsPath;

    public static ZooKeeper getZk() {
      return zk;
    }

    public static void setZk(ZooKeeper zk) {
      Config.zk = zk;
    }

    public static String getBoxNode() {
      return boxNode;
    }

    public static void setBoxNode(String boxNode) {
      Config.boxNode = boxNode;
    }


    public static String getJobsNodePath() {
      return jobsNodePath;
    }

    public static void setJobsNodePath(String jobsNodePath) {
      Config.jobsNodePath = jobsNodePath;
    }

    public static String getLeaderNodePath() {
      return leaderNodePath;
    }

    public static void setLeaderNodePath(String leaderNodePath) throws KeeperException,
        InterruptedException {
      Config.leaderNodePath = leaderNodePath;
      if (leaderNodePath.equals(getMyPath())) {
        Config.setIsLeader(true);
      } else {
        logger.info("I'm not the leader. My leader is " + leaderNodePath);
        Config.setIsLeader(false);
      }
    }


    public static String getMyPath() {
      return myPath;
    }

    public static void setMyPath(String myPath) {
      Config.myPath = myPath;
    }


    public static String getElectionsPath() {
      return electionsPath;
    }

    public static void setElectionsPath(String electionsPath) {
      Config.electionsPath = electionsPath;
    }

    private static boolean nodeBusy = false;

    // pointer to the job runner used in this to kill it when job is done.
    private static JobRunner runner = null;



    public static JobRunner getRunner() {
      return runner;
    }

    public static void setRunner(JobRunner runner) {
      Config.runner = runner;
    }

    public static boolean isNodeBusy() {
      return nodeBusy;
    }


    public static void setNodeBusy(boolean nodeBusy) {
      Config.nodeBusy = nodeBusy;
    }


    public static boolean getIsLeader() {
      return isLeader;
    }

    public static void setIsLeader(boolean isLeader) throws KeeperException, InterruptedException {
      if (isLeader) {
        logger.info("I just became leader! Yay! :)");
      }
      Config.isLeader = isLeader;
    }

    public static String getHostname() {
      return hostname;
    }

    private static void setHostname(String hostname) {
      Config.hostname = hostname;
    }
  }

  private static Logger logger = LoggerFactory.getLogger(ZooKeeperService.class);
  private static ObjectMapper jsonMapper = new ObjectMapper();

  @Inject
  Leader leader;

  private ZooKeeper connectToZooKeeper() throws IOException {
    logger.info("ZooKeeper Service Started! Welcome to the Zoo! :)");
    // try connecting to ZooKeeper
    String zkServer = SpringPropertiesUtil.getProperty("zk.server");
    String zkPort = SpringPropertiesUtil.getProperty("zk.port");
    logger.info("ZooKeeper Server: " + zkServer + " port: " + zkPort);
    return new ZooKeeper(zkServer + ":" + zkPort, 3000, this);
  }

  private void createBoxNode() throws KeeperException, InterruptedException {
    List<String> children = Config.getZk().getChildren(SEP, false);
    logger.info("Initial znodes in ZooKeeper root: " + children.toString());
    String boxId = SpringPropertiesUtil.getProperty("box.id");
    Config.setBoxNode(createNodeIfNotExists(SEP + boxId, new byte[0], Config.getZk(),
        Ids.OPEN_ACL_UNSAFE, CreateMode.PERSISTENT));
    logger.info("Should now have box id in znodes in ZooKeeper root: "
        + Config.getZk().getChildren(SEP, false));
  }

  private String registerOnElectionsNode() throws UnknownHostException, KeeperException,
      InterruptedException {
    // now there is a node for the box. lets create the elections node if it's not there
    Config.setElectionsPath(Config.getBoxNode() + SEP + ELECTIONS_NODE_STRING);
    createNodeIfNotExists(Config.getElectionsPath(), new byte[0], Config.getZk(),
        Ids.OPEN_ACL_UNSAFE, CreateMode.PERSISTENT);

    InetAddress addr = InetAddress.getLocalHost();
    Config.setHostname(addr.getHostName());

    String nodePath =
        Config.getZk().create(
            Config.getBoxNode() + SEP + ELECTIONS_NODE_STRING + SEP + NODE_PREFIX_STRING,
            Config.getHostname().getBytes(), Ids.OPEN_ACL_UNSAFE, CreateMode.EPHEMERAL_SEQUENTIAL);

    logger.info("Registered under the elections node. My node path is " + nodePath);
    return nodePath;
  }

  /**
   * {@link PostConstruct} method that is run right after the bean is initialized (only once)
   * 
   * @throws IOException
   * @throws KeeperException
   * @throws InterruptedException
   */
  @PostConstruct
  public synchronized void init() throws IOException, KeeperException, InterruptedException {

    Config.setZk(connectToZooKeeper());

    createBoxNode();

    Config.setJobsNodePath(Config.getBoxNode() + SEP + JOBS_NODE_STRING);

    Config.setMyPath(registerOnElectionsNode());

    Config.setLeaderNodePath(getNextLeader());

    // watch elections node
    Config.getZk().getChildren(Config.getBoxNode() + SEP + ELECTIONS_NODE_STRING, this);

    // remove all jobs from the jobs node...
    // TODO (amir) Add feature to be able to recover jobs that have been running after leaders
    // death. For now job queue is gone after leader dies for sanity reasons.

    if (Config.getIsLeader()) {
      try {
        logger
            .info("I'm the first leader... I will clean up the jobs node. This should be replaced"
                + " in future by a correct way of continuing executing from the queue from"
                + " previous leader.");
        deleteNodeRecursive(Config.getJobsNodePath());
      } catch (KeeperException ex) {
        if (ex.code().equals(Code.NONODE)) {
          // good. starting from a clean tree.
        }
      }
    }

    // also watch on jobs node or create it if it doesn't exist
    createNodeIfNotExists(Config.getJobsNodePath(), new byte[0], Config.getZk(),
        Ids.OPEN_ACL_UNSAFE, CreateMode.PERSISTENT);

    // get nodes and start watching!
    logger.info("Final leader: " + Config.getLeaderNodePath());
    watchNodeAndChildren(Config.getJobsNodePath());
  }


  public String createNodeIfNotExists(String path, byte[] data, ZooKeeper zk, List<ACL> acl,
      CreateMode createMode) throws KeeperException, InterruptedException {
    try {
      zk.create(path, data, acl, createMode);
      logger.debug("Created znode: " + path);
    } catch (KeeperException ex) {
      if (ex.code().equals(Code.NODEEXISTS)) {
        logger.debug("Did not create znode " + path + " because it already exists");
      }
    }
    return path;
  }

  public void processReelection() throws KeeperException, InterruptedException{
    logger.info("Election node has changed!");
    try {
      if (Config.getZk().exists(Config.getLeaderNodePath(), false) == null) {
        logger.warn("Previous leader was " + Config.getLeaderNodePath());
        logger
            .warn("My leader is dead! Will try to find a new leader or I will be the next leader if no one is left ;)");
        Config.setLeaderNodePath(getNextLeader());

        // watch elections node
        Config.getZk().getChildren(Config.getBoxNode() + SEP + ELECTIONS_NODE_STRING, this);
      } else {
        logger.info("Old leader leader (" + Config.getLeaderNodePath() + ") still alive.");
      }
    } catch (Exception e) {
      logger.error("Error in leader election!", e);
    }
    
    // keep watching elections node's children
    Config.getZk().getChildren(Config.getElectionsPath(), this);
  }
  
  public void processJobsNode() throws JsonParseException, JsonMappingException, KeeperException, InterruptedException, IOException{

    logger.info("Jobs node has changed!");
    // decide what to do
    // get the smallest job that needs a runner
    if (Config.getZk().getChildren(Config.getJobsNodePath(), false).size() > 0
        && !Config.isNodeBusy()) {
      logger.debug("Job's node path is: " + Config.getJobsNodePath());
      popAwaitingJobNode();
    }
    // re-register the box_id/jobs node watch
    watchNodeAndChildren(Config.getJobsNodePath());
  }
  
  public void processJobInitialization(WatchedEvent event, String lastJob) throws KeeperException, InterruptedException, NumberFormatException, IOException{
    String lastJobPath = Config.getJobsNodePath() + SEP + lastJob;
    logger.debug("CURRENT RUNNER NODE CHANGED: " + lastJobPath + SEP + RUNNERS_NODE_STRING);
    // TODO: Change this so jobs would have a min/max runner option
    JobRequest currentJobRequest = getJobRequest(lastJob);
    String runnersNode = Config.getJobsNodePath() + SEP + lastJob + SEP + RUNNERS_NODE_STRING;
    List<String> runners = Config.getZk().getChildren(runnersNode, false);
    logger.debug("Current job: " + lastJob + " started: " + currentJobRequest.getStarted()
        + " runners " + runners);

    if (runners.size() == Config.getZk().getChildren(Config.getElectionsPath(), false).size()
        && !currentJobRequest.getStarted()) {

      logger.info("There are enough runners for job " + lastJob
          + " therefore starting the run!");

      GenomeRegionSplitter splitter =
          new GenomeRegionSplitter(new File(
              currentJobRequest.getProperty(JobParameters.FASTA_FILE)), runners.size());
      GenomeBucket[] regions = splitter.getRegions();
      Map<GenomeBucket, String> genomeBucketToRunnerNodeMapping =
          new HashMap<GenomeBucket, String>();


      // create the genome bucket mapping
      for (int i = 0; i < runners.size(); i++) {
        RunnerNode runnerNode = getRunnerNode(runners.get(i), lastJob);
        genomeBucketToRunnerNodeMapping.put(regions[i], runnerNode.getId());
      }


      // mark job request as started
      currentJobRequest.setStarted(true);
      setJobRequest(currentJobRequest);

      // enough runners have registered. lets tell them to start!
      for (String runner : runners) {
        RunnerNode runnerNode = getRunnerNode(runner, lastJob);
        runnerNode.setStarted(true);
        runnerNode.setGenomeBucketToRunnerNodeMapping(genomeBucketToRunnerNodeMapping);

        // if it's a streaming job, wait for when all nodes are listening...
        if (currentJobRequest.isStreaming()) {
          watchNodeData(lastJobPath + SEP + RUNNERS_NODE_STRING + SEP + runner);
        }
        setRunnerNode(runnerNode, currentJobRequest.getId());
      }
    }

    // re-register watch to keep checking for job being done.
    watchNodeAndChildren(event.getPath());

  
  }
  
  public void processNodeChildrenChangedEvent(WatchedEvent event) throws KeeperException,
      InterruptedException, JsonParseException, JsonMappingException, IOException {
    if (event.getType() == EventType.NodeChildrenChanged) {
      
      // elections node changed
      if (event.getPath().equals(Config.getBoxNode() + SEP + ELECTIONS_NODE_STRING)) {
        processReelection();
      }

      // jobs node changed
      else if (event.getPath().equals(Config.getJobsNodePath())) {
        processJobsNode();
      }


      // need this for the next condition
      String lastJob =
          getNextSmallestNode(Config.getZk().getChildren(Config.getJobsNodePath(), false));

      // master checking if all nodes have registered or job is done
      if (Config.getIsLeader() && event.getPath().endsWith(SEP + RUNNERS_NODE_STRING)
          && lastJob != null) {
        processJobInitialization(event, lastJob);
      }

      processJobCompletion(event, lastJob);


    }

  }

  /**
   * @param event
   * @param lastJob
   * @throws InterruptedException
   */
  public void processJobCompletion(WatchedEvent event, String lastJob) throws InterruptedException {
    try {
      String lastJobPath = Config.getJobsNodePath() + SEP + lastJob;

      if (Config.getIsLeader()
          && event.getPath().equals(lastJobPath + SEP + DONE_RUNNERS_NODE_STRING)) {

        if (Config.getZk().getChildren(lastJobPath + SEP + DONE_RUNNERS_NODE_STRING, this).size() == Config
            .getZk().getChildren(Config.getElectionsPath(), false).size()
            && Config.getZk().getChildren(lastJobPath + SEP + RUNNERS_NODE_STRING, false).size() == 0) {
          // job seems to be be done
          logger.info("Job " + lastJobPath + " seems to be done! Deleting the job znode.");
          try {
            deleteNodeRecursive(lastJobPath);
          } catch (Exception e) {
            logger.error("Was not able to delete finished job node.", e);
          }
        } else {
          watchNodeAndChildren(lastJobPath + SEP + DONE_RUNNERS_NODE_STRING);
        }
      }
    } catch (KeeperException ex) {
      if (ex.code().equals(Code.NONODE)) {
        logger.debug("Job node has already been deleted...");
      }
    }
  }

  public void processNodeDataChanged(WatchedEvent event) throws Exception {

    if (event.getType() == EventType.NodeDataChanged) {
      String currentJob = getCurrentJobPath();
      JobRequest request = getJobRequest(parseNodeName(currentJob));


      if (event.getPath().equals(
          Config.getJobsNodePath() + SEP
              + getNextSmallestNode(Config.getZk().getChildren(Config.getJobsNodePath(), false))
              + SEP + RUNNERS_NODE_STRING + SEP + parseNodeName(Config.getMyPath()))) {

        // this is the runner node getting notified about data changed in it's node
        RunnerNode myNode =
            getRunnerNode(parseNodeName(Config.getMyPath()), getNextSmallestNode(Config.getZk()
                .getChildren(Config.getJobsNodePath(), false)));
        JobRequest jobRequest =
            getJobRequest(getNextSmallestNode(Config.getZk().getChildren(Config.getJobsNodePath(),
                false)));
        JobRunner runner;
        runner = new JobRunner(jobRequest, myNode);
        runner.start();
        Config.setRunner(runner);
      }



      // job has been KILLED!
      if (event.getPath().equals(currentJob) && request.isStopped()) {

        // this is the runner node getting notified about it's job being stopped!
        // kill the runner thread!

        Config.getRunner().interrupt();

        if (Config.getIsLeader()) {
          leader.getInputStreamer().interrupt();
        }
        Config.setNodeBusy(false);
        logger.info("Job " + request.getId() + " has been intrupted.");
      }



      if (Config.getIsLeader() && event.getPath().contains(JOBS_NODE_STRING)
          && event.getPath().contains(currentJob + SEP + RUNNERS_NODE_STRING)) {

        if (request.isStreaming() && request.getStarted() && !request.isStreamingStarted()) {
          // check if all the nodes have started listening.
          List<String> runners =
              Config.getZk().getChildren(currentJob + SEP + RUNNERS_NODE_STRING, false);
          if (runners.size() == Config.getZk().getChildren(Config.getElectionsPath(), false).size()) {
            for (String runner : runners) {
              try {
                RunnerNode runnerNode = getRunnerNode(runner, request.getId());
                runnerNode = getRunnerNode(runner, request.getId());
                runnerNode = getRunnerNode(runner, request.getId());
                if (!runnerNode.isListening()) {
                  logger
                      .debug("Not all the nodes have started listening... will keep watching the runners."
                          + event.getPath());
                  // keep watching
                  byte[] data = Config.getZk().getData(event.getPath(), false, null);
                  logger.debug("Value was: " + new String(data));

                  watchNodeData(event.getPath());
                  return;
                }
              } catch (KeeperException ex) {
                if (ex.code().equals(Code.NONODE)) {
                  // this watch went off because the node's data has been changed to done.
                  // just ignore (this is because I'm using NodeDataChanged trigger for both
                  // when nodes mark their runner as done and also when they start listening...
                  // I did this to keep the tree as simple as possible.
                  // this was probably not the best idea...

                  return;
                }
              }
            }


            request.setStreamingStarted(true);
            setJobRequest(request);

            // they are all LISTENING! :)
            logger.info("All runners are listening for job " + request.getId()
                + " will start streaming input.");


            logger.info("Will keep watching " + currentJob + " data in case the job gets killed.");
            // keep watching the job in case it is going to be killed!
            watchNodeData(currentJob);

            leader.streamInput(request);

          } else {
            logger.debug("Waiting for all the runners to start listening... " + runners.size()
                + " runners have started listening.");


          }
        } else {
          logger.warn("Why on earth did this get run?!");
        }
      }
    }
  }



  public void process(WatchedEvent event) {
    try {
      logger.info("Received this event from a watch: " + event.toString());
      //check if this is a NodeChildrenChanged event
      processNodeChildrenChangedEvent(event);
      //check if this is a NodeDataChanged event
      processNodeDataChanged(event);
    } catch (Exception e) {
      //TODO: some cleanup is needed here... but I dont know what yet. I'm sure I'll figure out as we go on...
      logger.error("Error processing ZooKeeper event " + event, e);
    }
  }

  @PreDestroy
  public void cleanUp() throws InterruptedException, KeeperException {
    if (Config.getZk() != null) {
      // remove my node from ZooKeeper
      logger.info("Deleting ZooKeeper ephemeral node before restarting the server.");
      Config.getZk().delete(Config.getMyPath(), 0);

      // remove all jobs from the jobs node...
      // TODO(amir) Add feature to be able to recover jobs that have been running after leaders
      // death. For now job queue is gone after leader dies for sanity reasons.
      if (Config.getIsLeader()) {
        deleteNodeRecursive(Config.getJobsNodePath());
      }

    } else {
      throw new RuntimeException("ZooKeeper Server doesn't exist!");
    }
  }

  private String getNextLeader() throws KeeperException, InterruptedException {
    // check if im the leader
    // criteria: if there is no other leader, i'm the one, otherwise the node with smallest id is
    // the leader
    List<String> currentNodes = Config.getZk().getChildren(Config.getElectionsPath(), false);

    String nextLeader = null;

    if (currentNodes.size() == 0) {
      // THIS SHOULD NOT HAPPEN
      logger
          .error("Something went wrong... getNextLeader() should be called when there has been a new node created.");
    } else if (currentNodes.size() == 1 && Config.getMyPath().contains(currentNodes.get(0))) {
      nextLeader = Config.getMyPath();
    } else {
      logger.debug("Current nodes are: " + currentNodes);
      nextLeader = Config.getElectionsPath() + SEP + getNextSmallestNode(currentNodes);
      logger.debug("Leader path: " + nextLeader);
    }
    return nextLeader;
  }

  public String getNextSmallestNode(List<String> children) {
    if (children.size() > 0) {
      long lowest = parseSequenceNumber(children.get(0));
      String lowestPath = parseNodeName(children.get(0));

      for (String child : children) {
        long current = parseSequenceNumber(child);
        if (current < lowest) {
          lowest = current;
          lowestPath = child;
        }
      }
      return lowestPath;
    }
    return null;
  }

  public long parseSequenceNumber(String znode) {
    return Integer.parseInt(znode.substring(znode.lastIndexOf(".") + 1));
  }

  public static String parseNodeName(String znodePath) {
    return znodePath.substring(znodePath.lastIndexOf("/") + 1);
  }

  private void watchNodeAndChildren(String nodePath) throws KeeperException, InterruptedException {
    Config.getZk().getChildren(nodePath, this);
  }

  private void watchNodeData(String nodePath) throws KeeperException, InterruptedException {
    Config.getZk().getData(nodePath, this, null);
  }

  /**
   * @param request
   * @throws Exception
   */
  public void createJob(JobRequest request) throws Exception {
    if (Config.getIsLeader()) {
      byte[] jsonJobConfig = jsonMapper.writeValueAsBytes(request);
      String jobsZnode = Config.getBoxNode() + SEP + JOBS_NODE_STRING + SEP;

      String job =
          Config.getZk().create(jobsZnode + JOB_NODE_STRING, jsonJobConfig, Ids.OPEN_ACL_UNSAFE,
              CreateMode.PERSISTENT_SEQUENTIAL);
      request.setId(parseNodeName(job));
      setJobRequest(request);

      // create runners node
      createNodeIfNotExists(jobsZnode + request.getId() + SEP + RUNNERS_NODE_STRING, new byte[0],
          Config.getZk(), Ids.OPEN_ACL_UNSAFE, CreateMode.PERSISTENT);

      // create done runners node
      createNodeIfNotExists(jobsZnode + request.getId() + SEP + DONE_RUNNERS_NODE_STRING,
          new byte[0], Config.getZk(), Ids.OPEN_ACL_UNSAFE, CreateMode.PERSISTENT);

      // start watching the job
      watchNodeAndChildren(jobsZnode + request.getId());

      // start watching the runners under job
      watchNodeAndChildren(jobsZnode + request.getId() + SEP + RUNNERS_NODE_STRING);

      // start watching the runnersdone under job

      watchNodeAndChildren(jobsZnode + request.getId() + SEP + DONE_RUNNERS_NODE_STRING);
    } else {
      throw new Exception("A non-leader called leader only method!");
    }
  }

  public int getDummyIp() {
    return (int) (19139 + parseSequenceNumber(Config.getMyPath()));
  }

  /**
   * @return
   */
  public static String getNodeName() {
    return parseNodeName(Config.getMyPath());
  }

  /**
   * @param job
   * @throws KeeperException
   * @throws InterruptedException
   * @throws IOException
   * @throws JsonMappingException
   * @throws JsonParseException
   */
  public static void runnerDone(JobRequest job) throws KeeperException, JsonParseException,
      JsonMappingException, IOException, InterruptedException {
    Stat stat = new Stat();
    Config.getZk().getData(
        Config.getJobsNodePath() + SEP + job.getId() + SEP + RUNNERS_NODE_STRING + SEP
            + parseNodeName(Config.getMyPath()), false, stat);
    Config.getZk().delete(
        Config.getJobsNodePath() + SEP + job.getId() + SEP + RUNNERS_NODE_STRING + SEP
            + parseNodeName(Config.getMyPath()), stat.getVersion());
    Config.getZk().create(
        Config.getJobsNodePath() + SEP + job.getId() + SEP + DONE_RUNNERS_NODE_STRING + SEP
            + parseNodeName(Config.getMyPath()), new byte[0], Ids.OPEN_ACL_UNSAFE,
        CreateMode.EPHEMERAL);
    Config.setNodeBusy(false);
    Config.getRunner().interrupt();
  }


  public static JobRequest getJobRequest(String jobId) throws KeeperException,
      InterruptedException, JsonParseException, JsonMappingException, IOException {
    String jobPath = Config.getBoxNode() + SEP + JOBS_NODE_STRING + SEP + jobId;
    byte[] data = Config.getZk().getData(jobPath, false, null);
    logger.debug("Looking up job " + jobPath);
    return jsonMapper.readValue(data, JobRequest.class);
  }

  public static void setJobRequest(JobRequest jobRequest) throws KeeperException,
      InterruptedException, JsonParseException, JsonMappingException, IOException {
    byte[] data = jsonMapper.writeValueAsBytes(jobRequest);
    Stat nodeStat = new Stat();
    Config.getZk().getData(Config.getBoxNode() + SEP + JOBS_NODE_STRING + SEP + jobRequest.getId(),
        false, nodeStat);
    Config.getZk().setData(Config.getBoxNode() + SEP + JOBS_NODE_STRING + SEP + jobRequest.getId(),
        data, nodeStat.getVersion());
  }

  public static RunnerNode getRunnerNode(String runnerId, String jobId) throws KeeperException,
      InterruptedException, JsonParseException, JsonMappingException, IOException {
    String runnerPath =
        Config.getBoxNode() + SEP + JOBS_NODE_STRING + SEP + jobId + SEP + RUNNERS_NODE_STRING
            + SEP + runnerId;
    Stat stat = new Stat();

    byte[] data = Config.getZk().getData(runnerPath, false, stat);
    logger.debug("Looking up runner " + runnerPath);

    SimpleModule module = new SimpleModule("EnhancedDatesModule", new Version(0, 1, 0, "alpha"));


    module.addKeyDeserializer(GenomeBucket.class, new GenomeBucketDeserializer());
    // and the magic happens here when we register module with mapper:
    jsonMapper.registerModule(module);

    return jsonMapper.readValue(data, RunnerNode.class);
  }

  public static void setRunnerNode(RunnerNode runnerNode, String jobId) throws KeeperException,
      InterruptedException, JsonParseException, JsonMappingException, IOException {
    byte[] data = jsonMapper.writeValueAsBytes(runnerNode);
    Stat nodeStat = new Stat();
    String runnerPath =
        Config.getBoxNode() + SEP + JOBS_NODE_STRING + SEP + jobId + SEP + RUNNERS_NODE_STRING
            + SEP + runnerNode.getId();
    Config.getZk().getData(runnerPath, false, nodeStat);
    Config.getZk().setData(runnerPath, data, nodeStat.getVersion());
  }

  private synchronized String popAwaitingJobNode() throws KeeperException, InterruptedException,
      JsonParseException, JsonMappingException, IOException {

    logger.debug("Popping the job with smallest id.");
    List<String> jobs = Config.getZk().getChildren(Config.getJobsNodePath(), false);


    String job = null;
    while (jobs.size() > 0) {
      job = getNextSmallestNode(jobs);
      try {
        JobRequest jobRequest = getJobRequest(job);
        if (!jobRequest.getStarted()) break;
      } catch (KeeperException ex) {
        if (ex.code().equals(Code.NONODE)) {
          logger.debug("Job got deleted when I was trying to add myself. Trying next one.");
        } else {
          throw ex;
        }
      }

      jobs.remove(job);
    }

    logger.debug("Job " + job + " seems to need a runner. Trying to add myself.");
    String jobNode = Config.getJobsNodePath() + SEP + job;

    // adding myself as a runner
    // get job config

    JobRequest jobRequest =
        jsonMapper.readValue(Config.getZk().getData(jobNode, false, null), JobRequest.class);
    jobRequest.setId(parseNodeName(jobNode));

    RunnerNode runnerNode = new RunnerNode();
    runnerNode.setHost(Config.getHostname());
    // TODO: Change these to read form Spring properties file
    runnerNode.setPort(getDummyIp());
    runnerNode.setSorterPort(getDummyIp() + 200);

    runnerNode.setId(parseNodeName(Config.getMyPath()));

    // add runner under runners for leader's info
    createNodeIfNotExists(jobNode + SEP + RUNNERS_NODE_STRING, new byte[0], Config.getZk(),
        Ids.OPEN_ACL_UNSAFE, CreateMode.PERSISTENT);


    // start watching runners node for job
    watchNodeAndChildren(jobNode + SEP + RUNNERS_NODE_STRING);

    Config.getZk().create(
        jobNode + SEP + RUNNERS_NODE_STRING + SEP + parseNodeName(Config.getMyPath()),
        jsonMapper.writeValueAsBytes(runnerNode), Ids.OPEN_ACL_UNSAFE, CreateMode.EPHEMERAL);

    // start watching my node on runners node
    watchNodeData(jobNode + SEP + RUNNERS_NODE_STRING + SEP + parseNodeName(Config.getMyPath()));

    Config.setNodeBusy(true);
    return jobNode;
  }

  public static List<RunnerNode> getJobRunners(JobRequest request) throws KeeperException,
      InterruptedException, JsonParseException, JsonMappingException, IOException {
    List<String> runners =
        Config.getZk().getChildren(
            Config.getJobsNodePath() + SEP + request.getId() + SEP + RUNNERS_NODE_STRING, false);
    List<RunnerNode> runnerNodes = new ArrayList<RunnerNode>();
    for (String runner : runners) {
      runnerNodes.add(getRunnerNode(runner, request.getId()));
    }
    return runnerNodes;
  }

  private String getCurrentJobPath() throws KeeperException, InterruptedException {
    String lastJob =
        getNextSmallestNode(Config.getZk().getChildren(Config.getJobsNodePath(), false));
    String lastJobPath = Config.getJobsNodePath() + SEP + lastJob;
    return lastJobPath;
  }

  private static void deleteNodeRecursive(String path) throws KeeperException, InterruptedException {
    List<String> children = Config.getZk().getChildren(path, false);
    Stat stat = new Stat();
    for (String child : children) {
      deleteNodeRecursive(path + "/" + child);
    }
    logger.info("Deleting " + path);
    Config.getZk().getData(path, false, stat);
    Config.getZk().delete(path, stat.getVersion());
  }

  public static void killJob(String jobId) throws JsonParseException, JsonMappingException,
      KeeperException, InterruptedException, IOException {
    JobRequest job = ZooKeeperService.getJobRequest(jobId);
    job.setStopped(true);
    ZooKeeperService.setJobRequest(job);
  }

  public static void killJobSilently(String jobId) {
    try {
      killJob(jobId);
    } catch (Exception ex) {
      logger
          .error(
              "Failed to stop failed job, i.e. disaster... Server might be in unusable state before a restart.",
              ex);
    }
  }

  public static boolean getIsLeader() {
    return Config.getIsLeader();
  }
}
