// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.io.IOException;
import java.net.InetAddress;
import java.util.List;

import javax.annotation.PostConstruct;
import javax.annotation.PreDestroy;

import org.apache.zookeeper.CreateMode;
import org.apache.zookeeper.KeeperException;
import org.apache.zookeeper.WatchedEvent;
import org.apache.zookeeper.Watcher;
import org.apache.zookeeper.Watcher.Event.EventType;
import org.apache.zookeeper.ZooDefs.Ids;
import org.apache.zookeeper.ZooKeeper;
import org.apache.zookeeper.data.ACL;
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
  private static final String JOBS_NODE_STRING = "jobs";
  private static final String NODE_PREFIX_STRING = "node.";
  private static InetAddress addr;
  private static String hostname;
  private static ZooKeeper zk;
  private static String boxNode;
  private static String jobsNodePath;
  private String leaderNodePath;
  static Integer mutex;
  private String myPath;
  private static String electionsPath;
  private static ObjectMapper jsonMapper = new ObjectMapper();
  
  @PostConstruct
  public void init() throws IOException, KeeperException, InterruptedException {
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

    getNextLeader();

    // get nodes and start watching!
    logger.info("Final leader: " + leaderNodePath);


    // also watch on jobs node or create it if it doesn't exist
    jobsNodePath = boxNode + SEP + JOBS_NODE_STRING;
    createNodeIfNotExists(jobsNodePath, new byte[0], zk, Ids.OPEN_ACL_UNSAFE, CreateMode.PERSISTENT);
    watchNodeAndChildren(jobsNodePath);

  }


  public String createNodeIfNotExists(String path, byte[] data, ZooKeeper zk, List<ACL> acl,
      CreateMode createMode) throws KeeperException, InterruptedException {
    if (zk.exists(path, false) == null) {
      zk.create(path, data, acl, createMode);
      logger.debug("Created znode: " + path);
    } else {
      logger.debug("Did not create znode " + path + " because it already exists");
    }
    return path;
  }

  public static boolean getIsLeader() {
    return isLeader;
  }

  public static void setIsLeader(boolean isLeader) {
    logger.info("I just became leader! Yay! :)");
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
      if (event.getType() == EventType.NodeChildrenChanged
          && event.getPath().equals(boxNode + SEP + ELECTIONS_NODE_STRING)) {
        logger.info("Election node has changed!");
        try {
          logger.warn("Previous leader was " + leaderNodePath);
          if (zk.exists(leaderNodePath, false) == null) {
            logger
                .warn("My leader is DEAD! Will try to find a new leader or I will be the next leader if no one is left ;)");
            getNextLeader();
            logger.warn("New leader is the her majesty: " + leaderNodePath + " my  path is "
                + myPath);
          } else {
            logger.warn("Elections node changed");
          }
        } catch (KeeperException e) {
          e.printStackTrace();
        } catch (InterruptedException e) {
          e.printStackTrace();
        }
      } else if (event.getType() == EventType.NodeChildrenChanged
          && event.getPath().equals(jobsNodePath)) {
        logger.info("Jobs node has changed!");
        // decide what to do
        // get the smallest job
        if (zk.getChildren(jobsNodePath, false).size() > 0) {
          String jobNode =
              jobsNodePath + SEP + getNextSmallestNode(zk.getChildren(jobsNodePath, false));

          // get job config
          JobRequest jobRequest =
              jsonMapper.readValue(zk.getData(jobNode, false, null),
                  JobRequest.class);
          JobRunner runner = new JobRunner();
          runner.run(jobRequest);

          // create runners node
          createNodeIfNotExists(jobNode + SEP + RUNNERS_NODE_STRING, new byte[0], zk,
              Ids.OPEN_ACL_UNSAFE, CreateMode.PERSISTENT);

          RunnerNode runnerNode = new RunnerNode();
          runnerNode.setHost(hostname);
          runnerNode.setIp(getDummyIp());

          // add runner under runners for leader's info
          zk.create(jobNode + SEP + RUNNERS_NODE_STRING + SEP + parseNodeName(myPath),
              jsonMapper.writeValueAsBytes(runnerNode), Ids.OPEN_ACL_UNSAFE,
              CreateMode.EPHEMERAL_SEQUENTIAL);
        }
      }



      if (boxNode != null) zk.getChildren(boxNode + SEP + ELECTIONS_NODE_STRING, this);
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
      isLeader = true;
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
    long lowest = parseSequenceNumber(myPath);
    String lowestPath = parseNodeName(myPath);
    for (String child : children) {
      long current = parseSequenceNumber(child);
      if (current < lowest) {
        lowest = current;
        lowestPath = child;
      }
    }
    return lowestPath;
  }

  public long parseSequenceNumber(String znode) {
    return Integer.parseInt(znode.substring(znode.lastIndexOf(".") + 1));
  }

  public String parseNodeName(String znodePath) {
    return znodePath.substring(znodePath.lastIndexOf("/") + 1);
  }

  private void watchNodeAndChildren(String nodePath) throws KeeperException, InterruptedException {
    zk.getChildren(nodePath, this);
  }


  /**
   * @param request
   * @throws Exception
   */
  public synchronized static void createJob(JobRequest request) throws Exception {
    if (getIsLeader()) {
      byte[] jsonJobConfig = jsonMapper.writeValueAsBytes(request);
      String jobsZnode = boxNode + SEP + JOBS_NODE_STRING + SEP;
      String jobZnode =
          zk.create(jobsZnode + request.getId(), jsonJobConfig, Ids.OPEN_ACL_UNSAFE,
              CreateMode.PERSISTENT_SEQUENTIAL);
    } else {
      throw new Exception("A non-leader called leader only method!");
    }
  }
  
  public long getDummyIp (){
    return 19139 + parseSequenceNumber(myPath);
  }


}
