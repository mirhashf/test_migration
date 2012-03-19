'''
Created on Apr 17, 2011

Stanford CS227b: General Game Playing
Assignment 3: Minimax Gamer
Team Rebecca Black Mesa
'''

from util.statemachine import MachineState
from util.statemachine.implementation.prover import ProverStateMachine

from player.gamer.statemachine import StateMachineGamer
from player.gamer.statemachine.reflex.event import ReflexMoveSelectionEvent

from useful import *
from time import time

class MinMaxGamer(StateMachineGamer):

    def __init__(self):
        StateMachineGamer.__init__(self)
        self.gameCache = {}
        self.MAX_DEPTH = 20
        self.THIS_DEPTH = 10
        self.TIMEOUT_BUFFER = 1900
        self.MAX_POSITIONS = 50000
        self.BF_EST = -1
        self.countStates = 0

    def setClassVariables(self):
        self.players = []
        roles = self.getStateMachine().getRoles()
        for i, role in enumerate(roles):
            if role == self.getRole():
                self.myIdx = i
            self.players.append( (i, role) )

    def getName(self):
        return "MinMaxGamer (Python)"
        
    def stateMachineMetaGame(self, timeout):
        pass
            
    def stateMachineSelectMove(self, timeout):
        #moves = self.getStateMachine().getLegalMoves(self.getCurrentState(), self.getRole())
        self.timeout = timeout
        self.setClassVariables()
        state = self.getCurrentState()
        self.countStates = 0
        
        playerMoves = map(lambda (idx, role) : self.getStateMachine().getLegalMoves(state, role), self.players)
        if len(playerMoves[self.myIdx]) == 1:
            return playerMoves[self.myIdx][0]
            
        allMoves = cartesian_product(*playerMoves)        
        allMoves = list(allMoves)
        
        num_moves = len(playerMoves[self.myIdx])
        if self.BF_EST == -1:
            self.BS_EST = num_moves
        else:
            self.BS_EST = self.BS_EST + num_moves / 2
        depth = 1
        while (num_moves**(depth+1) < self.MAX_POSITIONS):
            depth = depth + 1
        self.THIS_DEPTH = depth
        
        try:
            children = [(i, self.minimax(self.getStateMachine().getNextState(state, move), min, 0)) for i, move in enumerate(allMoves)]
            selectionIdx = max(children, key=lambda x:x[1])[0]
            selection = allMoves[selectionIdx][self.myIdx]

        except TimeoutException:
            print "TIMEOUT"            
            selection = playerMoves[self.myIdx][0]

        print self.countStates

        self.notifyObservers(ReflexMoveSelectionEvent(playerMoves[self.myIdx], selection, 1))
        return selection
        
    def stateMachineStop(self):
        pass
        
    def getInitialStateMachine(self):
        return ProverStateMachine()
                        
                
    def minimax(self, state, func, depth):    
        self.countStates += 1
        if self.timeout - time()*1000 < self.TIMEOUT_BUFFER:
            raise TimeoutException                
                
        hashCode = state.hashCode()
        if hashCode in self.gameCache:
            return self.gameCache[hashCode]
    
        if self.getStateMachine().isTerminal(state):       
            self.gameCache[hashCode] = self.getStateMachine().getGoal(state, self.getRole())
            return self.gameCache[hashCode]
                        
        if depth >= self.MAX_DEPTH:
            return 25        
                           
        playerMoves = map(lambda (idx, role) : self.getStateMachine().getLegalMoves(state, role), self.players)
        allMoves = cartesian_product(*playerMoves)
        children = [self.getStateMachine().getNextState(state, move) for move in allMoves]

        next_level_func = (min if func == max else max)
        values = map(lambda child: self.minimax(child, next_level_func, depth+1), children)
        self.gameCache[hashCode] = func(values)
        return self.gameCache[hashCode]
        