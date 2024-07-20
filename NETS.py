"""
Author: Alexander Söderhäll
Date: 2024-02-24
Description:    NETS Outlier detection, translated from Java GitHub project: https://github.com/kaist-dmlab/NETS
                They are the original creators of NETS.
"""

import math


"""
@class          Tuple is the data class in NETS, effectively a Datapoint class with some extra values.       

@param          id: Unique ID of this tuple.
@param          slideID: Which slide in arrived at.
@param          value: Datapoint object.
"""
class Tuple:
    def __init__(self, id: int, slideID: int, value: float):
        self.id = id
        self.slideID = slideID
        self.value = value
        self.unSafeOutNeighbors = {}
        self.lastNNSlideID = -1
        self.safeness = False

        self.nn = 0
        self.nnSafeOut = 0
        self.nnUnSafeOut = 0
        self.nnIn = 0
        self.fullDimCellIdx = 0
        self.subDimCellIdx = 0


    """
    @function       Returns how many neighbors that the tuple has.
    """
    def getNN(self):
        nn = self.nnIn + self.nnSafeOut + self.nnUnSafeOut
        return nn


    """
    @function       Removes outliers that exit the window.

    @param          itr: Current timestep of SUMO.
    @param          nS: Current slide amount.
    """
    def removeOutdatedNNUnsafeOut(self, itr: int, nS: int):
        for key in self.unSafeOutNeighbors:
            if key <= itr-nS:
                self.nnUnSafeOut -= self.unSafeOutNeighbors[key]


    """
    @function       Remove all unsafeoutneighbors. Never used.
    """
    def truncate(self):
        self.unSafeOutNeighbors.clear()
    

    """
    @function       Simple print for a Tuple.
    """
    def printIt(self):
        print("ID:", self.id, "Datapoint", self.value)


"""
@class          A cell is a datastructure that NETS use to perform outlier detecion. A cell has Tuple (datapoints) and a lenght of at most R.

@param          cellIndex: Index of this cell.
@param          cellCenter: Center of this cell.
"""
class Cell:
    def __init__(self, cellIndex: int, cellCenter: float):
        self.cellCenter = cellCenter
        self.tuples = set()
        self.cellIdx = cellIndex


    """
    @function       Returns the amount of tuples in the cell.
    """
    def getNumTuples(self):
        return len(self.tuples)


    """
    @function       Adds a tuple object to the cell.

    @param          tuple: A tuple object.
    """
    def addTuple(self, tuple: Tuple):
        self.tuples.add(tuple)


"""
@class          NETS is conducts distance-based outlier detection in ONE dimension only. To extend this implementation, please follow GitHub Java
                https://github.com/kaist-dmlab/NETS/blob/master for details on what should be expanded.

                
@param          dim: Dimenstion of data.
@param          subDim: subdimension of data.
@param          R: Radius which dictates if another Tuple is a neighbor.
@param          R: Radius which dictates if another Tuple is a neighbor.
@param          k: Minimum nighbors in order to be inlier.
@param          S: Slide value.
@param          W: Window size.
@param          nW: ??
@param          maxValues: max value for each dimension.
@param          minValues: min value for each dimension.
"""
class NETS:
    def __init__(self, dim: int, subDim: int, R: float, k: int, S: int, W: int, nW: int, maxValues: list, minValues: list):

        # Neighborhood range R.
        self.R = R

        # Minimum amount of neighbors k to be an inlier.
        self.k = k

        # Window slide value w.slide.
        self.S = S

        # Window size w.s.
        self.W = W

        # Never used here nor in Java implementation.
        self.nW = nW

        # Amount of slides that go into one window.
        self.nS = W/S

        self.subDim = subDim
        self.dim = dim
        self.subDimFlag = False if subDim == dim else True

        # We only implement NETS for one dimension, otherwise first is sqrt(sumDim)*2 and second sqrt(dim)*2
        self.neighCellIdxDist = math.sqrt(self.subDim)*2
        self.neighCellFullDimIdxDist = math.sqrt(self.dim)*2

        # In our case, speed limit of an edge.
        self.maxValues = maxValues

        # Minimum speed is 0 m/s.
        self.minValues = minValues


        minDimSize = float('inf')
        dimSize = [0]*self.dim

        for i in range(self.dim):
            dimSize[i] = maxValues[i] - minValues[i]
            if dimSize[i] < minDimSize:
                minDimSize = dimSize[i]

        dimWeightsSum = 0
        dimWeights = [0]*self.dim
        for i in range(self.dim):
            dimWeights[i] = 1
            dimWeightsSum += dimWeights[i]
        
        self.dimLength = [0]*self.dim
        gapCount = [0]*self.dim

        for i in range(self.dim):
            self.dimLength[i] = math.sqrt(self.R**2*dimWeights[i]/dimWeightsSum)
            gapCount[i] = math.ceil(dimSize[i]/self.dimLength[i])
            dimSize[i] = gapCount[i]*self.dimLength[i]

        if self.subDimFlag:
            minSubDimSize = float('inf')
            subDimSize = [0]*self.subDim
            for i in range(self.subDim):
                subDimSize[i] = self.maxValues[i] - self.minValues[i]
                if subDimSize[i] < minSubDimSize:
                    minSubDimSize = subDimSize[i]

            subDimWeightsSum = 0
            subDimWeights = [0]*self.subDim
            for i in range(self.subDim):
                subDimWeights[i] = 1
                subDimWeightsSum += subDimWeights[i]
            
            self.subDimLength = [0]*self.subDim
            subDimgapCount = [0]*self.subDim
            for i in range(self.subDim):
                self.subDimLength[i] = math.sqrt(self.R**2*subDimWeights[i]/subDimWeightsSum)
                subDimgapCount[i] = math.ceil(subDimSize[i]/self.subDimLength[i])
                subDimSize[i] = subDimgapCount[i]*self.subDimLength[i]

        


        self.slideIn = {}
        
        self.windowCnt = {}
        self.slides = []
        self.slideOut = {}
        self.idxDecoder = {}
        self.idxEncoder = {}
        self.fullDimCellWindowCnt = {}
        self.fullDimCellSlidesCnt = []
        self.fullDimCellSlideOutCnt = {}
        self.candidateCellsTupleCnt = 0

        self.fullDimCellSlideInCnt = {}
        self.slideDelta = {}

        self.influencedCells = set()

        self.outliers = set()


    """
    @function       Calculates the amount of outliers.

    @returns        Returns the amount of outliers.
    """
    def showOutliers(self):
        return len(self.outliers)


    """
    @function       Returns all outliers in a list, their values (data point objects) only.

    @returns        List of Datapoint objects which are outliers.
    """
    def returnOutliers(self):
        dataPointlist = []
        for elem in self.outliers:
            dataPointlist.append(elem.value)
        return dataPointlist


    """
    @function       Given a list of newly arrived tuples, inserts them into the correct cellindex.
                    If this is a first tuple in an index, initialise a new Cell.

    @param          slideTuples: A list of newly arrived Tuples.
    """
    def indexingSlide(self, slideTuples: list):
        self.slideIn = {}
        self.fullDimCellSlideInCnt = {}

        for t in slideTuples:
            fullDimCellIdx = []
            subDimCellIdx = []
            for j in range(self.dim):
                dimIdx = t.value[j] - self.minValues[j] / self.dimLength[j]
                fullDimCellIdx.append(dimIdx)
            if self.subDimFlag:
                for j in range(self.subDim):
                    dimIdx = t.value[j] - self.minValues[j] / self.subDimLength[j]
                    subDimCellIdx.append(dimIdx)
            else:
                subDimCellIdx = fullDimCellIdx

            t.fullDimCellIdx = fullDimCellIdx
            t.subDimCellIdx = fullDimCellIdx
            if tuple(fullDimCellIdx) not in self.idxEncoder:
                id = len(self.idxEncoder)
                self.idxEncoder[tuple(fullDimCellIdx)] = id
                self.idxDecoder[id] = fullDimCellIdx
            
            if tuple(subDimCellIdx) not in self.idxEncoder:
                id = len(self.idxEncoder)
                self.idxEncoder[tuple(subDimCellIdx)] = id
                self.idxDecoder[id] = subDimCellIdx

            if self.idxEncoder[tuple(fullDimCellIdx)] not in self.slideIn:
                cellCenter = [0.0]*self.subDim
                #print("cellCenter", cellCenter)
                if self.subDimFlag:
                    for j in range(self.subDim):
                        cellCenter[j] = self.minValues[j] + subDimCellIdx[j]*self.subDimLength[j]+self.subDimLength[j]/2
                else:
                    for j in range(self.subDim):
                        cellCenter[j] = self.minValues[j] + fullDimCellIdx[j]*self.dimLength[j]+self.dimLength[j]/2
                self.slideIn[self.idxEncoder[tuple(fullDimCellIdx)]] = Cell(cellIndex=fullDimCellIdx, cellCenter=cellCenter)
                """ print("fullDim",fullDimCellIdx)
                print("tuple value", tuple.value) """

            self.slideIn[self.idxEncoder[tuple(fullDimCellIdx)]].addTuple(t)
            
            #print("slideIn", self.slideIn)

            if self.idxEncoder[tuple(fullDimCellIdx)] not in self.fullDimCellSlideInCnt:
                self.fullDimCellSlideInCnt[self.idxEncoder[tuple(fullDimCellIdx)]] = 0
            
            self.fullDimCellSlideInCnt[self.idxEncoder[tuple(fullDimCellIdx)]] += 1
        
        self.slides.append(self.slideIn)
        #print(len(self.slides))
        self.fullDimCellSlidesCnt.append(self.fullDimCellSlideInCnt)
        #print(self.idxEncoder)
        #print("decoder", self.idxDecoder, "\n")


    """
    @function       Removes old slide and adds data from new slide. Also checks if this leads to a change in the amount of tuples in a cell index.

    @param          slideTuples: A list of Tuples (Tuple.value is a Datapoint object).
    @param          itr: The current timestep of the SUMO simulation. 
    """
    def calcNetChange(self, slideTuples: list, itr: int):
        #print(slideTuples)

        """ for elem in slideTuples:
            print("Datapoint uid, timestamp", elem.value.uid, elem.value.timestamp) """
        self.indexingSlide(slideTuples)
        """ if self.edge == self.attackEdge:
            print(len(self.slides), itr) """
        
        # Each time this function is called, we've reached a full window, so we want to pop the prior slide.
        if itr > self.W:
            """ if self.edge == self.attackEdge:
                print("got here") """
            self.slideOut = self.slides.pop(0)
            self.fullDimCellSlideOutCnt = self.fullDimCellSlidesCnt.pop(0)

        #Old
        """ if itr > self.nS-1:
            self.slideOut = self.slides.pop(0)
            self.fullDimCellSlideOutCnt = self.fullDimCellSlidesCnt.pop(0) """
        self.slideDelta = {}


        for key in self.slideIn:
            if key not in self.windowCnt:
                self.windowCnt[key] = 0
                self.idxDecoder[key] = self.slideIn[key].cellIdx
                #print("slideIn key", key)
                #print("cellIdx", self.slideIn[key].cellIdx)
            
            self.windowCnt[key] = self.windowCnt[key] + self.slideIn[key].getNumTuples()
            self.slideDelta[key] = self.slideIn[key].getNumTuples()

        for key in self.slideOut:
            self.windowCnt[key] = self.windowCnt[key] - self.slideOut[key].getNumTuples()
            if self.windowCnt[key] < 1:
                del self.windowCnt[key]
            
            if key in self.slideDelta:
                self.slideDelta[key] = self.slideDelta[key] - self.slideOut[key].getNumTuples()
            else:
                self.slideDelta[key] = self.slideOut[key].getNumTuples()*-1
            
        for key in self.fullDimCellSlideInCnt:
            if key not in self.fullDimCellWindowCnt:
                self.fullDimCellWindowCnt[key] = 0
            self.fullDimCellWindowCnt[key] +=  self.fullDimCellSlideInCnt[key]
        
        for key in self.fullDimCellSlideOutCnt:
            self.fullDimCellWindowCnt[key] -= self.fullDimCellSlideOutCnt[key]
            if self.fullDimCellWindowCnt[key] < 1:
                del self.fullDimCellWindowCnt[key]
        #print(self.slideDelta)
        #print("windowCnt", self.windowCnt)
        #print(self.slides)


    """
    @function       Given an update, will retrieve all cell indexes which have had a change in the amount of data in their cell.
    """
    def getInfCellIndices(self):
        self.influencedCells = set()
        for cellIdxWin in self.windowCnt:
            if self.windowCnt[cellIdxWin] > self.k:
                continue
            """ if self.edge == self.attackEdge:
                print("cellIDxWin", cellIdxWin) """
            for cellIdxSld in self.slideDelta:
                if self.neighboringSet(self.idxDecoder[cellIdxWin], self.idxDecoder[cellIdxSld]):
                    """ if self.edge == self.attackEdge:
                        print(cellIdxWin, cellIdxSld) """
                    if cellIdxWin not in self.influencedCells:
                        self.influencedCells.add(cellIdxWin)
                    break
        """ if self.edge == self.attackEdge:
            print(self.influencedCells) """


    """
    @function       Given the index of an influenced cell, retrieve all cells which neighbor that cell and
                    could be influenced by the influenced cell.

    @param          cellIdxInf: The index of the influenced cell.

    @returns        Returns a list of cells that may be influenced by the newly updated cell with index @cellIdxInf.
    """
    def getSortedCandidateCellIndices(self, cellIdxInf: int):
        candidateCellIndices = []

        candidateCellIndicesMap = {}
        for cellIdxWin in self.windowCnt:
            dist = self.neighboringSetDist(self.idxDecoder[cellIdxInf], self.idxDecoder[cellIdxWin])
            #print("Dist", dist)
            if cellIdxInf != cellIdxWin and dist < self.neighCellIdxDist:
                if dist not in candidateCellIndicesMap:
                    candidateCellIndicesMap[dist] = [cellIdxWin]
        #print("candidateCellIndicesMap",candidateCellIndicesMap)
        keys = sorted(candidateCellIndicesMap.keys())
        for key in keys:
            candidateCellIndices.extend(candidateCellIndicesMap[key])
            for cellIdxWin in candidateCellIndicesMap[key]:
                #print(cellIdxWin)
                #print("windowCnt", self.windowCnt[cellIdxWin])
                self.candidateCellsTupleCnt += self.windowCnt[cellIdxWin]
        
        return candidateCellIndices


    """
    @function       Since new data has entered, check the slide that is being removed and remove those tuples from necessary 
                    data structures. Then find outliers.

    @param          itr: The current timestep iteration (time in the MCS system).
    """
    def findOutlier(self, itr: int):
        outliersToRemove = set()
        #print(self.slideOut)
        
        """ if(slideOut.containsKey(idxEncoder.get(outlier.subDimCellIdx)) && slideOut.get(idxEncoder.get(outlier.subDimCellIdx)).tuples.contains(outlier)) {  
				it.remove();
			}else if(fullDimCellWindowCnt.get(idxEncoder.get(outlier.fullDimCellIdx))>K){ 
				it.remove();
		} """

        # For each outlier, check if they are leaving (@self.slideOut) and if so remove them.
        for outlier in self.outliers:
            if(self.idxEncoder[tuple(outlier.subDimCellIdx)] in self.slideOut and outlier in self.slideOut[self.idxEncoder[outlier.subDimCellIdx]].tuples):
                #print("Got here 1")
                outliersToRemove.add(outlier)
            elif self.fullDimCellWindowCnt[self.idxEncoder[tuple(outlier.fullDimCellIdx)]] > self.k:
                #print("Got here 2")
                outliersToRemove.add(outlier)

        # Update the set, remove @outliersToRemove.
        self.outliers.difference_update(outliersToRemove)
        
        self.findOutlierNETS(itr)
        """ if self.edge == self.attackEdge:
            print("DATA\n")
            #print("slideIn",self.slideIn)
            print("windowCnt",self.windowCnt)
            #print("slides",self.slides)
            print("slideOut",self.slideOut)
            print("idxDecoder",self.idxDecoder)
            print("idxEncoder",self.idxEncoder)
            print("fullDimCellWindowCnt",self.fullDimCellWindowCnt)
            print("fullDimCellSlidesCnt",self.fullDimCellSlidesCnt)
            print("fullDimCellSlideOutCnt",self.fullDimCellSlideOutCnt)
            print("candidateCellsTupleCnt",self.candidateCellsTupleCnt)
            print("fullDimCellSlideInCnt",self.fullDimCellSlideInCnt)
            print("slideDelta",self.slideDelta)
            print("influencedCells",self.influencedCells)
            print("\nEND") """
        
        
    """
    @function       Finds outliers with the NETS implementation. Starts by getting influenced cells (those that either had data removed or inserted when new slide arrived).    
                    Then 

    @param          itr: The current timestep iteration (time in the MCS system).
    """
    def findOutlierNETS(self, itr: int):
        self.getInfCellIndices()
        #print("Amount of slides:", len(self.windowCnt))
        exitFully = False
        """ if self.edge == self.attackEdge:
            print("Inf cells",self.influencedCells) """
        while(True):
            # Java code uses GOTO statements, this does the same.
            if exitFully:
                break

            # For each influenced cell, update them and check if they are now outliers (they do not have enough tuples.)
            for infCellIdx in self.influencedCells:
                """ if self.edge == self.attackEdge:
                    print("awdwad", self.influencedCells, infCellIdx) """
                self.candidateCellsTupleCnt = 0
                candCellIndices = self.getSortedCandidateCellIndices(infCellIdx)
                """ if self.edge == self.attackEdge:
                    print(candCellIndices) """
                """ if self.edge == self.attackEdge:
                    print("windowCnt[infCellIdx]", self.windowCnt[infCellIdx])
                    print("self.candidateCellsTupleCnt",self.candidateCellsTupleCnt) """
                self.candidateCellsTupleCnt += self.windowCnt[infCellIdx]
                """ if self.edge == self.attackEdge:    
                    print("candidateCellsTupleCnt", self.candidateCellsTupleCnt) """
                if self.candidateCellsTupleCnt < self.k + 1:
                    """ if self.edge == self.attackEdge:
                        print(infCellIdx) """
                    for slide in self.slides:
                        if infCellIdx not in slide:
                            continue
                        self.outliers.update(slide[infCellIdx].tuples)
                        #print(self.outliers)
                    #continue
                    continue
                
                # For the influenced cells (those that had data augmented from new slide entering)
                # update properties of each tuple, so they are now inliers.
                candOutlierTuples = set()
                for slide in self.slides:
                    if infCellIdx not in slide:
                        continue
                    for t in slide[infCellIdx].tuples:
                        if t.safeness:
                            continue
                        t.nnIn = self.fullDimCellWindowCnt[self.idxEncoder[tuple(t.fullDimCellIdx)]] - 1
                        t.removeOutdatedNNUnsafeOut(itr, self.nS)
                        #print(tuple.getNN())
                        if t.getNN() < self.k:
                            candOutlierTuples.add(t)
                        elif t in self.outliers:
                            #print(tuple.value)
                            self.outliers.remove(t)

                # For all cadidate outliers (those that have less than k neighbors)
                # check neighboring cells if those tuples could help.
                while(True):
                    TupleLoop = False
                    #print("aaaaaa")
                    #print(candOutlierTuples)
                    for tCand in candOutlierTuples:
                        if TupleLoop:
                            break
                        currentSlideID = itr+1
                        #print("bbbbb")
                        while(True):
                            #print("Im here!")
                            for currentSlide in reversed(self.slides):
                                if TupleLoop:
                                    break
                                currentSlideID -= 1
                                if currentSlideID in tCand.unSafeOutNeighbors:
                                    break
                                else:
                                    tCand.unSafeOutNeighbors[currentSlideID] = 0
                                
                                while(True):
                                    for otherCellIdx in candCellIndices:
                                        #print("otherCellIdx", otherCellIdx)
                                        #print("candCellIndices", candCellIndices)
                                        if TupleLoop:
                                            break
                                        if otherCellIdx not in currentSlide or not self.neighboringTupleSet(tCand.value, currentSlide[otherCellIdx].cellCenter, 1.5*self.R):
                                            break

                                        otherTuples = set()

                                        otherTuples = currentSlide[otherCellIdx].tuples
                                    
                                        for tOther in otherTuples:
                                            if self.neighboringTuple(tCand, tOther, self.R):
                                                if tCand.slideID <= tOther.slideID:
                                                    tCand.nnSafeOut += 1
                                                else:
                                                    tCand.nnUnSafeOut += 1
                                                    tCand.unSafeOutNeighbors[currentSlideID] += 1
                                                if tCand.nnSafeOut >= self.k:
                                                    if(tCand in self.outliers):
                                                        #print(tCand.value)
                                                        self.outliers.remove(tCand)
                                                        tCand.safeness = True
                                                    TupleLoop = True
                                                    break
                                    if tCand.getNN() >= self.k:
                                        if tCand in self.outliers:
                                            #print(tCand.value)
                                            self.outliers.remove(tCand)
                                        TupleLoop = True
                                        break
                                    break
                        
                            self.outliers.add(tCand)
                            #print("Got here")
                            break
                    break
            exitFully = True
            break


    """
    @function       Determines the L2-norm between to tuples values (data point speed), detemines if they are neighbors.

    @param          t1: First tuple.
    @param          t2: Second tuple.
    @param          threshold: R.

    @returns        Returns True if @ss is below the threshold, false otherwise.
    """
    def neighboringTuple(self, t1, t2, threshold):
        ss = 0.0

        # Since we're taking the L2-norm, multiply R by itself so it's in the same "scale".
        threshold = threshold * threshold
        ss += (t1.value - t2.value)**2
        if ss > threshold:
            return False
        return True


    """
    @function       Takes a tuple @v1 in one cell and a cell centre @v2 from another cell and checks if they are within
                    a distance of 1.5*R.

    @param          v1: Tuple value 1.
    @param          v2: Cell center of a tuple in another cell.
    @param          threshold: 1.5*R

    @returns        Returns True if @ss is below the threshold, false otherwise.
    """
    def neighboringTupleSet(self, v1, v2, threshold):
        ss = 0.0
        threshold = threshold**2
        for k in range(len(v2)):
            ss += (v1[k]-v2[k])**2
            if ss > threshold:
                return False
        return True


    """
    @function       Compares the relative cell index distance.

    @param          c1: One cell index.
    @param          c2: Another cell index.

    @returns        Returns inf distance if @ss is above the threshold, the l2 norm otherwise.
    """
    def neighboringSetDist(self, c1, c2):
        ss = 0.0
        cellIdxDist = self.neighCellFullDimIdxDist if len(c1) == self.dim else self.neighCellIdxDist
        threshold = cellIdxDist**2
        for k in range(len(c1)):
            ss += (c1[k]-c2[k])**2
            if ss >= threshold:
                return float('inf')
        return math.sqrt(ss)
    

    """
    @function       Compares the relative cell index distance and checks if they are within distance of dim (2R^2).

    @param          c1: One cell index.
    @param          c2: Another cell index.

    @returns        Returns False is @ss is greater than some threshold, True otherwise.
    """
    def neighboringSet(self, c1, c2):
        ss = 0.0
        cellIdxDist = self.neighCellFullDimIdxDist if len(c1) == self.dim else self.neighCellIdxDist
        threshold = cellIdxDist**2
        for k in range(len(c1)):
            ss += (c1[k]-c2[k])**2
            if ss >= threshold:
                return False
        return True