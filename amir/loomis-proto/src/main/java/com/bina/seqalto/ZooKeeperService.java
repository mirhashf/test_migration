// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.io.IOException;
import java.net.InetAddress;
import java.util.List;

import javax.annotation.PostConstruct;
import javax.annotation.PreDestroy;

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
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Service;

/**
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 * 
 */
@Service
public class ZooKeeperService implements Watcher {
  private static Logger logger = LoggerFactory.getLogger(ZooKeeperService.class);
  private static boolean isLeader;
  private static final String SEP = "/";
  private static final String ELECTIONS_NODE_STRING = "elections";
  private static final String RUNNERS_NODE_STRING = "runners";
  private static final String DONE_RUNNERS_NODE_STRING = "done_runners";
  private static final String JOBS_NODE_STRING = "jobs";
  private static final String JOB_NODE_STRING = "job.";
  private static final String NODE_PREFIX_STRING = "node.";
  private static InetAddress addr;
  private static String hostname;
  private static ZooKeeper zk;
  private static String boxNode;
  private static String jobsNodePath;
  private String leaderNodePath;
  static Integer mutex;
  private static String myPath;
  private static String electionsPath;
  private static ObjectMapper jsonMapper = new ObjectMapper();
  private static boolean nodeBusy = false;
  private static String previousLeader;

  @PostConstruct
  public synchronized void init() throws IOException, KeeperException, InterruptedException {
    logger.info("ZooKeeper Service Started! Welcome to the Zoo! :)");
    // try connecting to ZooKeeper
    String zkServer = SpringPropertiesUtil.getProperty("zk.server");
    String zkPort = SpringPropertiesUtil.getProperty("zk.port");
    logger.info("ZooKeeper Server: " + zkServer + " port: " + zkPort);
    zk = new ZooKeeper(zkServer + ":" + zkPort, 3000, this);

    mutex = new Integer(-1);

    List<String> children = zk.getChildren(SEP, false);
    logger.info("Initial znodes in ZooKeeper root: " + children.toString());
    String boxId = SpringPropertiesUtil.getProperty("box.id");
    boxNode =
        createNodeIfNotExists(SEP + boxId, new byte[0], zk, Ids.OPEN_ACL_UNSAFE,
            CreateMode.PERSISTENT);
    logger.info("Should now have box in znodes in ZooKeeper root: " + children.toString());

    // now there is a node for the box. lets create the elections node if it's not there
    electionsPath = boxNode + SEP + ELECTIONS_NODE_STRING;
    createNodeIfNotExists(electionsPath, new byte[0], zk, Ids.OPEN_ACL_UNSAFE,
        CreateMode.PERSISTENT);

    addr = InetAddress.getLocalHost();
    setHostname(addr.getHostName());


    myPath =
        zk.create(boxNode + SEP + ELECTIONS_NODE_STRING + SEP + NODE_PREFIX_STRING, getHostname()
            .getBytes(), Ids.OPEN_ACL_UNSAFE, CreateMode.EPHEMERAL_SEQUENTIAL);

    logger.info("My path is " + myPath);
    jobsNodePath = boxNode + SEP + JOBS_NODE_STRING;

    // also watch on jobs node or create it if it doesn't exist
    createNodeIfNotExists(jobsNodePath, new byte[0], zk, Ids.OPEN_ACL_UNSAFE, CreateMode.PERSISTENT);
    getNextLeader();

    // get nodes and start watching!
    logger.info("Final leader: " + leaderNodePath);

    watchNodeAndChildren(jobsNodePath);

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

  public static boolean getIsLeader() {
    return isLeader;
  }

  public void setIsLeader(boolean isLeader) throws KeeperException, InterruptedException {
    if (isLeader) {
      logger.info("I just became leader! Yay! :)");
    }
    ZooKeeperService.isLeader = isLeader;
  }


  public static String getHostname() {
    return hostname;
  }


  private static void setHostname(String hostname) {
    ZooKeeperService.hostname = hostname;
  }


  public synchronized void process(WatchedEvent event) {
    try {
      logger.info("Received this event from a watch: " + event.toString());
      if (event.getType() == EventType.NodeChildrenChanged) {

        // elections node changed
        if (event.getPath().equals(boxNode + SEP + ELECTIONS_NODE_STRING)) {
          logger.info("Election node has changed!");
          try {
            if (zk.exists(leaderNodePath, false) == null) {
              logger.warn("Previous leader was " + leaderNodePath);
              previousLeader = leaderNodePath;
              logger
                  .warn("My leader is DEAD! Will try to find a new leader or I will be the next leader if no one is left ;)");
              getNextLeader();
              logger.warn("New leader is the her majesty: " + leaderNodePath + " my path is "
                  + myPath);
            } else {
              logger.info("Old leader leader (" + leaderNodePath + ") still alive.");
            }
          } catch (KeeperException e) {
            e.printStackTrace();
          } catch (InterruptedException e) {
            e.printStackTrace();
          }
        }

        // jobs node changed
        else if (event.getPath().equals(jobsNodePath)) {
          logger.info("Jobs node has changed!");
          // decide what to do
          // get the smallest job that needs a runner
          if (zk.getChildren(jobsNodePath, false).size() > 0 && !isNodeBusy()) {
            logger.debug("Job's node path is: " + jobsNodePath);
            popAwaitingJobNode();

          }
          // re-register the box_id/jobs node watch
          watchNodeAndChildren(jobsNodePath);
        }

        String lastJob = getNextSmallestNode(zk.getChildren(jobsNodePath, false));
        String lastJobPath = jobsNodePath + SEP + lastJob;

        logger.debug("is Leader: " + getIsLeader());

        // master checking if all nodes have registered or job is done
        if (getIsLeader() && event.getPath().endsWith(SEP + RUNNERS_NODE_STRING)) {
          logger.debug("CURRENT RUNNER NODE CHANGED: " + lastJobPath + SEP + RUNNERS_NODE_STRING);
          // TODO: Change this so jobs would have a min/max runner option
          JobRequest currentJobRequest = getJobRequest(lastJob);
          String runnersNode = jobsNodePath + SEP + lastJob + SEP + RUNNERS_NODE_STRING;
          List<String> runners = zk.getChildren(runnersNode, false);
          logger.debug("Current job: " + lastJob + " started: " + currentJobRequest.getStarted()
              + " runners " + runners);

          if (runners.size() == zk.getChildren(electionsPath, false).size()
              && !currentJobRequest.getStarted()) {

            logger.info("There are enough runners for job " + lastJob
                + " therefore starting the run!");

            // mark job request as started
            currentJobRequest.setStarted(true);
            setJobRequest(currentJobRequest);

            // enough runners have registered. lets tell them to start!
            for (String runner : runners) {
              RunnerNode runnerNode = getRunnerNode(runner, lastJob);
              runnerNode.setStarted(true);
              setRunnerNode(runnerNode, lastJob);
            }

          }
          // re-register watch to keep checking for job being done.
          watchNodeAndChildren(event.getPath());

        }

        try {
          if (getIsLeader() && event.getPath().equals(lastJobPath + SEP + DONE_RUNNERS_NODE_STRING)) {

            if (zk.getChildren(lastJobPath + SEP + DONE_RUNNERS_NODE_STRING, this).size() == zk
                .getChildren(electionsPath, false).size()
                && zk.getChildren(lastJobPath + SEP + RUNNERS_NODE_STRING, false).size() == 0) {
              // job seems to be be done
              logger.info("Job " + lastJobPath + " seems to be done! Deleting the job znode.");
              try {
                deleteNodeRecursive(lastJobPath);
              } catch (Exception e) {
                e.printStackTrace();
                logger.error("Was not able to delete finished job node.");
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


      if (event.getType() == EventType.NodeDataChanged) {
        if (event.getPath().equals(
            jobsNodePath + SEP + getNextSmallestNode(zk.getChildren(jobsNodePath, false)) + SEP
                + RUNNERS_NODE_STRING + SEP + parseNodeName(myPath))) {

          // this is the runner node getting notified about data changed in it's node
          RunnerNode myNode =
              getRunnerNode(parseNodeName(myPath),
                  getNextSmallestNode(zk.getChildren(jobsNodePath, false)));
          JobRunner jobRunner = new JobRunner();
          JobRequest jobRequest =
              getJobRequest(getNextSmallestNode(zk.getChildren(jobsNodePath, false)));
          jobRunner.run(jobRequest, myNode);

        }


        if (boxNode != null) zk.getChildren(boxNode + SEP + ELECTIONS_NODE_STRING, this);
      }
    } catch (KeeperException e) {
      e.printStackTrace();
    } catch (InterruptedException e) {
      e.printStackTrace();
    } catch (JsonParseException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (JsonMappingException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }


  @PreDestroy
  public void cleanUp() throws InterruptedException, KeeperException {
    if (zk != null) {
      // remove my node from ZooKeeper
      logger.info("Deleting ZooKeeper ephemeral node before restarting the server.");
      zk.delete(myPath, 0);
    } else {
      throw new RuntimeException("ZooKeeper Server doesn't exist!");
    }
  }

  private void getNextLeader() throws KeeperException, InterruptedException {
    // check if im the leader
    // criteria: if there is no other leader, i'm the one, otherwise the node with smallest id is
    // the leader
    List<String> currentNodes = zk.getChildren(boxNode + SEP + ELECTIONS_NODE_STRING, false);

    if (currentNodes.size() == 0) {
      // THIS SHOULD NOT HAPPEN
      logger
          .error("Something went wrong... getNextLeader() should be called when there has been a new node created.");
    } else if (currentNodes.size() == 1 && myPath.contains(currentNodes.get(0))) {
      setIsLeader(true);
      leaderNodePath = myPath;
    } else {
      logger.debug("Current nodes are: " + currentNodes);
      leaderNodePath = electionsPath + SEP + getNextSmallestNode(currentNodes);
      logger.debug("Next leader path: " + leaderNodePath);
      if (leaderNodePath.equals(myPath)) {
        setIsLeader(true);
      } else {
        logger
            .info("I'm not the leader. I'm just an insignificant follower.One day I'll be a leader...");
        setIsLeader(false);
      }
    }
    zk.getChildren(boxNode + SEP + ELECTIONS_NODE_STRING, this);
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
    zk.getChildren(nodePath, this);
  }

  private void watchNodeData(String nodePath) throws KeeperException, InterruptedException {
    zk.getData(nodePath, this, null);
  }

  /**
   * @param request
   * @throws Exception
   */
  public void createJob(JobRequest request) throws Exception {
    if (getIsLeader()) {
      byte[] jsonJobConfig = jsonMapper.writeValueAsBytes(request);
      String jobsZnode = boxNode + SEP + JOBS_NODE_STRING + SEP;

      String job =
          zk.create(jobsZnode + JOB_NODE_STRING, jsonJobConfig, Ids.OPEN_ACL_UNSAFE,
              CreateMode.PERSISTENT_SEQUENTIAL);
      request.setId(parseNodeName(job));
      setJobRequest(request);

      // create runners node
      createNodeIfNotExists(jobsZnode + request.getId() + SEP + RUNNERS_NODE_STRING, new byte[0],
          zk, Ids.OPEN_ACL_UNSAFE, CreateMode.PERSISTENT);

      // create done runners node
      createNodeIfNotExists(jobsZnode + request.getId() + SEP + DONE_RUNNERS_NODE_STRING,
          new byte[0], zk, Ids.OPEN_ACL_UNSAFE, CreateMode.PERSISTENT);

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

  public long getDummyIp() {
    return 19139 + parseSequenceNumber(myPath);
  }

  /**
   * @return
   */
  public static String getNodeName() {
    return parseNodeName(myPath);
  }

  /**
   * @param job
   * @throws KeeperException
   * @throws InterruptedException
   * @throws IOException
   * @throws JsonMappingException
   * @throws JsonParseException
   */
  public static void runnerDone(JobRequest job) throws InterruptedException, KeeperException,
      JsonParseException, JsonMappingException, IOException {
    Stat stat = new Stat();
    zk.getData(jobsNodePath + SEP + job.getId() + SEP + RUNNERS_NODE_STRING + SEP
        + parseNodeName(myPath), false, stat);
    zk.delete(jobsNodePath + SEP + job.getId() + SEP + RUNNERS_NODE_STRING + SEP
        + parseNodeName(myPath), stat.getVersion());
    zk.create(jobsNodePath + SEP + job.getId() + SEP + DONE_RUNNERS_NODE_STRING + SEP
        + parseNodeName(myPath), new byte[0], Ids.OPEN_ACL_UNSAFE, CreateMode.EPHEMERAL);
    setNodeBusy(false);
  }


  public static boolean isNodeBusy() {
    return nodeBusy;
  }


  public static void setNodeBusy(boolean nodeBusy) {
    ZooKeeperService.nodeBusy = nodeBusy;
  }

  public static JobRequest getJobRequest(String jobId) throws KeeperException,
      InterruptedException, JsonParseException, JsonMappingException, IOException {
    String jobPath = boxNode + SEP + JOBS_NODE_STRING + SEP + jobId;
    byte[] data = zk.getData(jobPath, false, null);
    logger.debug("Looking up job " + jobPath);
    return jsonMapper.readValue(data, JobRequest.class);
  }

  public static void setJobRequest(JobRequest jobRequest) throws KeeperException,
      InterruptedException, JsonParseException, JsonMappingException, IOException {
    byte[] data = jsonMapper.writeValueAsBytes(jobRequest);
    Stat nodeStat = new Stat();
    zk.getData(boxNode + SEP + JOBS_NODE_STRING + SEP + jobRequest.getId(), false, nodeStat);
    zk.setData(boxNode + SEP + JOBS_NODE_STRING + SEP + jobRequest.getId(), data,
        nodeStat.getVersion());
  }

  public static RunnerNode getRunnerNode(String runnerId, String jobId) throws KeeperException,
      InterruptedException, JsonParseException, JsonMappingException, IOException {
    String runnerPath =
        boxNode + SEP + JOBS_NODE_STRING + SEP + jobId + SEP + RUNNERS_NODE_STRING + SEP + runnerId;
    byte[] data = zk.getData(runnerPath, false, null);
    logger.debug("Looking up runner " + runnerPath);
    return jsonMapper.readValue(data, RunnerNode.class);
  }

  public static void setRunnerNode(RunnerNode runnerNode, String jobId) throws KeeperException,
      InterruptedException, JsonParseException, JsonMappingException, IOException {
    byte[] data = jsonMapper.writeValueAsBytes(runnerNode);
    Stat nodeStat = new Stat();
    String runnerPath =
        boxNode + SEP + JOBS_NODE_STRING + SEP + jobId + SEP + RUNNERS_NODE_STRING + SEP
            + runnerNode.getId();
    zk.getData(runnerPath, false, nodeStat);
    zk.setData(runnerPath, data, nodeStat.getVersion());
  }

  private synchronized String popAwaitingJobNode() throws KeeperException, InterruptedException,
      JsonParseException, JsonMappingException, IOException {

    logger.debug("Popping the job with smallest id.");
    List<String> jobs = zk.getChildren(jobsNodePath, false);


    String job = null;
    while (jobs.size() > 0) {
      job = getNextSmallestNode(jobs);
      JobRequest jobRequest = getJobRequest(job);
      if (!jobRequest.getStarted()) break;
      jobs.remove(job);
    }

    logger.debug("Job " + job + " seems to need a runner. Trying to add myself.");
    String jobNode = jobsNodePath + SEP + job;

    // adding myself as a runner
    // get job config

    JobRequest jobRequest =
        jsonMapper.readValue(zk.getData(jobNode, false, null), JobRequest.class);
    jobRequest.setId(parseNodeName(jobNode));

    RunnerNode runnerNode = new RunnerNode();
    runnerNode.setHost(hostname);
    runnerNode.setIp(getDummyIp());
    runnerNode.setId(parseNodeName(myPath));

    // add runner under runners for leader's info
    createNodeIfNotExists(jobNode + SEP + RUNNERS_NODE_STRING, new byte[0], zk,
        Ids.OPEN_ACL_UNSAFE, CreateMode.PERSISTENT);
    

    // start watching runners node for job
    watchNodeAndChildren(jobNode + SEP + RUNNERS_NODE_STRING);
    
    zk.create(jobNode + SEP + RUNNERS_NODE_STRING + SEP + parseNodeName(myPath),
        jsonMapper.writeValueAsBytes(runnerNode), Ids.OPEN_ACL_UNSAFE, CreateMode.EPHEMERAL);

    // start watching my node on runners node
    watchNodeData(jobNode + SEP + RUNNERS_NODE_STRING + SEP + parseNodeName(myPath));
    

    
    setNodeBusy(true);
    return jobNode;
  }

  private static void deleteNodeRecursive(String path) throws KeeperException, InterruptedException {
    List<String> children = zk.getChildren(path, false);
    Stat stat = new Stat();
    for (String child : children) {
      deleteNodeRecursive(path + "/" + child);

    }
    logger.info("Deleting " + path);
    zk.getData(path, false, stat);
    zk.delete(path, stat.getVersion());
  }
}
