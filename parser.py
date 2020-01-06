# python3 parser.py -i ../HG002_PacBio_GRCh38.bam -r ../GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  -o test.txt
# python3 parser.py -i ../HG002.10kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam -r ../GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  -o test.txt 1> teststdout.vcf 2> teststderr.vcf

import argparse
import sys
import math
import pysam
from graphviz import Digraph
from pyfasta import Fasta
from Bio import pairwise2
from itertools import combinations
from scipy import stats
from collections import defaultdict
import copy

import logging
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

def findStartFromCigar(read, start, ref=None, seq=None, tIndex=None):
    """
    Find the start index of a sequence given the position in reference and cigar
    :param read: a pysam aligned segment
    :param start: position in reference
    :param ref: saved ref position for skipping forward in the cigar
    :param seq: saved sequence position for skipping forward in the cigar
    :param tIndex: index to skip to in the cigar
    """
    if ref:
        refPos = ref
    else:
        refPos = read.reference_start

    if seq:
        seqPos = seq
    else:
        seqPos = 0

    if tIndex is None:
        tIndex = 0

    for t in range(tIndex, len(read.cigartuples)):
        # Move reference genome position based on cigar tuple
        if read.cigartuples[t][0] in [0, 2, 3, 7, 8]:  # Only advance if match, deletion, skip, equal, or mismatch
            if (refPos + read.cigartuples[t][1]) > start:
                savedRef = refPos
                savedSeq = seqPos
                if read.cigartuples[t][0] in [0, 7, 8]:
                    # save the ref and seq here so when finding the end can start from there
                    while refPos < start:
                        refPos += 1
                        seqPos += 1
                return [read, refPos, t, seqPos-1, savedRef, savedSeq]
            refPos += read.cigartuples[t][1]
        if read.cigartuples[t][0] in [0, 1, 4, 7, 8]:  # Only advance if match, softclip, insertion, equal, or mismatch
            seqPos += read.cigartuples[t][1]
        if read.cigartuples[t][0] in [5, 6]:
            print(read.cigartuples[t][1], read.cigartuples[t][0])
    return [read, refPos, t, seqPos - 1, refPos, seqPos]

def findEndFromCigar(read, refPos, tIndex, seqPos, end):
    """
    Find the ending index in the read corresponding to the end position based on the cigar
    :param read: a pysam aligned segment
    :param refPos: saved ref position for skipping forward in the cigar
    :param tIndex: index to skip to in the cigar
    :param seqPos: saved sequence position for skipping forward in the cigar
    :param end: end reference position
    """
    for t in range(tIndex, len(read.cigartuples)):
        # Move reference genome position based on cigar tuple
        if read.cigartuples[t][0] in [0, 2, 3, 7, 8]:  # Only advance if match, deletion, skip, equal, or mismatch
            if (refPos + read.cigartuples[t][1]) > end:
                if read.cigartuples[t][0] in [0, 7, 8]:
                    while refPos < end:
                        refPos += 1
                        seqPos += 1
                return [read, refPos, t, seqPos]
            refPos += read.cigartuples[t][1]
        if read.cigartuples[t][0] in [0, 1, 4, 7, 8]:  # Only advance if match, insertion, equal, or mismatch
            seqPos += read.cigartuples[t][1]
        if read.cigartuples[t][0] in [5, 6]:
            print(read.cigartuples[t][1], read.cigartuples[t][0])
    # if end past end of sequence
    return [read, refPos, t, seqPos]



def parseBam(inputBam, inputRef, outputFile):
    """
    Goes through bam files and calls variants.
    :param inputBam: the bam file to call variants on
    :param inputRef: the fasta reference file
    :param outputFile: currently not used
    """
    samfile = pysam.AlignmentFile(inputBam, "rb")
    outFile = open(outputFile, "w+")

    chr = "chr20"
    # startIndex = 0
    # endIndex = 64444167
    # startIndex = 130270
    # endIndex = startIndex + 30
    # startIndex = 1019800 # no reads are seen
    # endIndex = startIndex + 30
    # startIndex = 2092640 #short repeat
    # endIndex = startIndex + 30
    # startIndex = 2145400 #long repeat
    # endIndex = startIndex + 30
    # startIndex = 2145400
    # endIndex = startIndex + 10000
    # startIndex = 2105150
    # endIndex = startIndex + 30


    startIndex = 2000000
    endIndex = 3000000
    k = 8
    windowSize = 30
    windowOverlap = 5

    variantDict = {}
    allVariantsDict = {}

    if endIndex - startIndex < windowSize:
        endIndex = startIndex + windowSize
        print("endIndex needs to be larger than or equal to startIndex + windowSize. Setting endIndex to: ", endIndex)

    print("##fileformat=VCFv4.3")
    print("##reference=GCA_000001405.15")
    print('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">')
    print('##INFO=<ID=FA,Number=1,Type=Integer,Description="Fraction of samples following GT genotype">')
    print('##INFO=<ID=FF,Number=1,Type=Integer,Description="Fraction of reads in the forward direction">')
    print('##INFO=<ID=AS,Number=1,Type=Integer,Description="Sum of alignment score for one genotype subtracted by the other">')
    print('##INFO=<ID=NW,Number=1,Type=Integer,Description="Number of windows variant was found in">')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    if endIndex- startIndex < 10000:
        fetchSize = endIndex - startIndex
    else:
        fetchSize = 10000
    for largeSeg in range(startIndex, endIndex - windowSize+1, fetchSize):
        largeSegStart = largeSeg
        largeSegEnd = largeSeg + fetchSize
        cursors = []
        reads = samfile.fetch('chr20', largeSegStart, largeSegEnd)
        savedRead = next(reads, None)


        if savedRead != None:

            for interval in range(largeSegStart, largeSegEnd - windowSize + 1, windowOverlap):
                index1 = interval
                index2 = interval + windowSize
                myGraph = Graph(chr + "_" + str(index1) + "_" + str(index2), chr, index1, index2, k)

                # get the ref sequence
                f = Fasta(inputRef)
                ref = f["chr20  AC:CM000682.2  gi:568336004  LN:64444167  rl:Chromosome  M5:b18e6c531b0bd70e949a7fc20859cb01  AS:GRCh38"][index1-1:index2]

                changedk = k
                i = index1
                refNodes = deBruijn(myGraph, changedk, ref, firstPos=index1, ref=True)
                # while refNodes == False and changedk <= 15:
                #     changedk += 1
                #     myGraph = Graph(chr + "_" + str(index1) + "_" + str(index2), chr, index1, index2, changedk)
                #     refNodes = deBruijn(myGraph, changedk, ref, ref=True)
                # if refNodes == False:
                #     refNodes = deBruijn(myGraph, changedk, ref, ref=True, final=True)

                myGraph.ref = ref

                # check if need to add any new cursors
                while savedRead:
                    read = savedRead
                    refPos = read.reference_start
                    if refPos >= index2:
                        break
                    cursors.append(findStartFromCigar(read, index1))
                    savedRead = next(reads, None)

                removeIndexes = []

                # update cursors and add seq to graph
                savedWindowedReads = []
                for i in range(0, len(cursors)):
                    read, refPosStart, cigarIndex, seqPosStart, savedRef, savedSeq = cursors[i]
                    cursors[i] = findStartFromCigar(read, index1, ref=savedRef, seq=savedSeq, tIndex=cigarIndex)
                    read, refPosStart, cigarIndex, seqPosStart, savedRef, savedSeq = cursors[i]
                    # if refPosStart >= index1 and refPosStart <= index2:
                    # get the appropiate sequence
                    r, refPosEnd, t, seqPosEnd = findEndFromCigar(read, savedRef, cigarIndex, savedSeq, index2)
                    # print(seqPosEnd, index2)
                    # if refPosEnd <= index2 and refPosStart <= index1:
                    if r.reference_end > index1 and r.reference_start < index2:
                        windowedRead = read.query_alignment_sequence[seqPosStart:seqPosEnd]
                        if len(windowedRead) and len(windowedRead) < windowSize *2:
                            print(windowedRead)
                            savedWindowedReads.append(windowedRead)
                            deBruijn(myGraph, changedk, windowedRead, reverse=read.is_reverse)

                    else:
                        # remove from memory
                        removeIndexes.append(i)
                for index in removeIndexes[::-1]:
                    del cursors[index]

                startingNodes = myGraph.pruneGraph(2)
                myGraph.multiplicity(savedWindowedReads)

                myGraph.printGraph()

                # find paths
                if len(startingNodes) == 0:
                    startingNodes.append(refNodes[0])
                myGraph.performSearch(startingNodes)

                # make sure to add the ref seq if doesn't exist.
                refPath = -1
                for i in range(0, len(myGraph.discovered)):
                    if myGraph.discovered[i] == refNodes:
                        refPath = i
                if refPath == -1:
                    myGraph.discovered.append(refNodes)
                    refPath = len(myGraph.discovered) - 1


                myGraph.refPath = refPath

                alignments = []
                variants4eachPath = []
                # which path do you belong to the best?
                for path in myGraph.discovered:
                    newAlignments = []

                    for read in savedWindowedReads:
                        newAlignments.append(max(i[2] for i in pairwise2.align.globalms(mergeNodes(path), read,  2, -1, -1.5, -.1)))
                        # print(pairwise2.align.globalms(mergeNodes(path), read, 2, -1, -.5, -.1))
                    alignments.append(newAlignments)

                    variant4path = {}
                    maxPathAlign = []
                    maxScore = -math.inf
                    for i in pairwise2.align.globalms(ref, mergeNodes(path), 2, -1, -1.5, -.1):
                        if maxScore < i[2]:
                            maxPathAlign = i
                            maxScore = i[2]

                    # TODO only calculate variants for the best two paths
                    for REF, ALT, POS, PATH in findDifference(maxPathAlign[0], maxPathAlign[1], index1, path):
                        if len(PATH) > 1:
                            variant4path = myGraph.storeVariants(POS, REF, ALT, PATH, variant4path)
                    variants4eachPath.append(variant4path)



                # instead of finding the best path, we want to find the reference path.
                refSumCount = sum(alignments[refPath])
                outFile.write("position: "+ str(index1) + "\t numPaths" + str(len(myGraph.discovered))+"\n")
                if len(myGraph.discovered) > 100: #TODO get rid of this arbitrary parameter
                    outFile.write("0\n")
                elif len(myGraph.discovered) > 1:
                    for i in range(0, len(myGraph.discovered)):
                        if i == refPath:
                            outFile.write("\t"+mergeNodes(myGraph.discovered[i])+"ref"+"\n")
                        else:
                            outFile.write("\t"+mergeNodes(myGraph.discovered[i])+"\n")
                    best2Comb = []
                    best2AvgAlignment = -1
                    best2AlignmentForEachRead = []
                    best2Paths = []
                    comb = combinations([i for i in range(0, len(myGraph.discovered))], 2)

                    # for each combination
                    for c in comb:
                        pathForEachRead = []
                        alignmentForEachRead = []

                        # for each read
                        for i in range(0, len(savedWindowedReads)):
                            bestPathForRead = -1
                            bestAlignmentScoreForRead = -1

                            # for each path in combination
                            for j in c:
                                if bestAlignmentScoreForRead < alignments[j][i]:
                                    bestPathForRead = j
                                    bestAlignmentScoreForRead = alignments[j][i]
                            pathForEachRead.append(bestPathForRead)
                            alignmentForEachRead.append(bestAlignmentScoreForRead)
                        if len(alignmentForEachRead):
                            averageAlign = sum(alignmentForEachRead) / len(alignmentForEachRead)
                        else:
                            averageAlign = 0
                        if averageAlign > best2AvgAlignment:
                            best2Comb = c
                            best2AvgAlignment = averageAlign
                            best2Paths = pathForEachRead
                            best2AlignmentForEachRead = alignmentForEachRead

                    comb1 = []
                    comb2 = []
                    for i in range(0, len(best2AlignmentForEachRead)):
                        if best2Paths[i] == best2Comb[0]:
                            comb1.append(best2AlignmentForEachRead[i])
                        else:
                            comb2.append(best2AlignmentForEachRead[i])



                    bestPathForEachRead = [comb1, comb2]

                    if refPath in best2Comb:
                        outFile.write("2\n")
                    else:
                        # TODO will not print if there is only one alt path with all the reads
                        outFile.write("3\n")

                    # find the variants and add to currentPositionVariants
                    variantDict = myGraph.findVariants(best2Comb, bestPathForEachRead, best2AlignmentForEachRead, refSumCount, variantDict, variants4eachPath)
                    variantDict = myGraph.printVariants(variantDict, index1)

                    allVariantsDict = myGraph.findVariants([i for i in range(0, len(myGraph.discovered))],bestPathForEachRead, best2AlignmentForEachRead,refSumCount, allVariantsDict, variants4eachPath)
                    allVariantsDict = myGraph.printVariants(allVariantsDict, index1, debug=True)
                else:
                    outFile.write("1\n")


    # print the rest of the variants that haven't been removed from the list
    # do this at the end or else ordering might be off
    myGraph.printVariants(variantDict)
    myGraph.printVariants(allVariantsDict, debug=True)
    outFile.write("\n")
    samfile.close()

def findClosestIndex(refPos, index, lessThan=True):
    """
    Search through the refPos list to find the closest position to index
    """
    lastRealRefPos = 0
    i = 0
    for pos in refPos:
        if pos is not None:
            if pos > index:
                if lessThan:
                    if i == 0 or lastRealRefPos == 0:
                        return 0
                    return lastRealRefPos-1
                else:
                    return lastRealRefPos
            lastRealRefPos = i
        i += 1
    return -1


class Node:
    """
    Create nodes in a graph and store them in a dictionary.

    Attributes:
        self.out - edges going out
        self.ins - edges going in
        self.name - name of node
    """
    def __init__(self, name, ref=False):
        """Store the node attributes and add node to dictionary of nodes."""
        self.out = {}
        self.ins = {}
        self.name = name
        self.ref = ref
        self.refPos = []
        self.visited = False
        self.inDegree = 0
        self.readPos = []
        self.refMultiplicity = 0
        self.multiplicity = []

    def addIn(self, newIn, reverse, ref):
        """Add a new in to node."""
        addCount = 1
        if ref:
            addCount = 0
        if newIn in self.ins:
            self.ins[newIn][reverse] += addCount
        else:
            self.ins[newIn] = [0, 0]
            self.ins[newIn][reverse] += addCount

    def addOut(self, newOut, reverse, ref):
        """Add a new out and in to node, if haven't seen out before, add to dictionary."""
        addCount = 1
        if ref:
            addCount = 0
        if newOut in self.out:
            self.out[newOut][reverse] += addCount
        else:
            self.out[newOut] = [0, 0]
            self.out[newOut][reverse] += addCount

    def removeOut(self, toRemoveOut, reverse, ref):
        addCount = -1
        if ref:
            addCount = 0
        if toRemoveOut in self.out:
            self.out[toRemoveOut][reverse] += addCount
        if self.out[toRemoveOut] == [0, 0]:
            del toRemoveOut

    def removeIn(self, toRemoveIn, reverse, ref):
        addCount = -1
        if ref:
            addCount = 0
        if toRemoveIn in self.out:
            self.ins[toRemoveIn][reverse] += addCount
        if self.ins[toRemoveIn] == [0, 0]:
            del toRemoveIn

    def addRefPos(self, pos):
        """
        Add a reference pos to a ref node
        """
        if pos not in self.refPos:
            self.refPos.append(pos)

    def addReadPos(self, pos):
        """
        Add a read pos to a node
        """
        self.readPos.append(pos)

class Variant:
    """
    Store a variant.
    """
    def __init__(self, CHROM, POS, ID, REF, ALT, MINCOUNT, FRACFORWARD):
        self.CHROM = CHROM
        self.POS = POS
        self.ID = ID
        self.REF = REF
        self.ALT = ALT
        self.NS = MINCOUNT
        self.FF = [FRACFORWARD]
        self.AS = [-1]
        self.FA = [-1]
        self.NW = 1


class Graph:
    """
    A collection of nodes.
    """
    def __init__(self, name, chr, pos1, pos2, k):
        self.name = name
        self.nodes = defaultdict(list)
        self.totalEdges = 0
        self.ref = ""
        self.chr = chr
        self.pos1 = pos1
        self.pos2 = pos2
        self.k = k
        self.discovered = []
        self.refPath = ""

    def addEdge(self, strRep1, strRep2, node1, node2, reverse, ref):
        """
        Add an edge from node1 to node2
        """
        if strRep1 not in self.nodes:
            self.nodes[strRep1].append(node1)
        elif ref:
            self.nodes[strRep1][0].addRefPos(node1.refPos[0])
        if strRep2 not in self.nodes:
            self.nodes[strRep2].append(node2)
        elif ref:
            self.nodes[strRep2][0].addRefPos(node2.refPos[0])

        if len(self.nodes[strRep1]) == 1 and len(self.nodes[strRep2]) == 1:
            self.nodes[strRep1][0].addOut(self.nodes[strRep2][0], reverse, ref)
            self.nodes[strRep2][0].addIn(self.nodes[strRep1][0], reverse, ref)
        # if check:
        #     if self.topologicalOrder() == False:
        #         # remove nodes
        #         # create a new Node - create a new Node class?
        #         # or add to Node if at correct position


    def refAddEdge(self, strRep1, strRep2, node1, node2, reverse, ref, pos):
        if node1 not in self.nodes[strRep1]:
            self.nodes[strRep1].append(node1)
        if node2 not in self.nodes[strRep2]:
            self.nodes[strRep2].append(node2)

        node1.addOut(node2, reverse, ref)
        node2.addIn(node1, reverse, ref)
        node1.addReadPos(pos)
        node2.addReadPos(pos+1)
        # else:
        #     check = True
        #     if strRep1 not in self.nodes:
        #         check = False
        #         if node1 not in self.nodes[strRep1]:
        #             self.nodes[strRep1].append(node1)
        #     if strRep2 not in self.nodes:
        #         check = False
        #         if node2 not in self.nodes[strRep2]:
        #             self.nodes[strRep2].append(node2)
        #
        #     if check == False:
        #         node1.addOut(node2, reverse, ref)
        #         node2.addIn(node1, reverse, ref)
        #         node1.addReadPos(pos)
        #         node2.addReadPos(pos + 1)
        #     else:
        #
        #
        #
        #         if len(self.nodes[strRep1]) == 1 and len(self.nodes[strRep2]) == 1:
        #             self.nodes[strRep1][0].addOut(self.nodes[strRep2][0], reverse, ref)
        #             self.nodes[strRep2][0].addIn(self.nodes[strRep1][0], reverse, ref)
        #             if self.topologicalOrder() == False:
        #                 self.nodes[strRep1][0].removeOut(self.nodes[strRep2][0], reverse, ref)
        #                 self.nodes[strRep2][0].removeIn(self.nodes[strRep1][0],reverse, ref)
        #                 self.nodes[strRep2].append(node2)
        #                 self.nodes[strRep1][0].addOut(node2, reverse, ref)
        #                 node2.addIn(self.nodes[strRep1][0], reverse, ref)
        #
        #                 self.nodes[strRep1][0].addReadPos(pos)
        #                 node2.addReadPos(pos + 1)
        #                 # if self.topologicalOrder() == False:
        #                 #     self.nodes[strRep1][0].removeOut(node2, reverse, ref)
        #                 #     node2.removeIn(self.nodes[strRep1][0], reverse, ref)
        #                 #     node1.addOut(self.nodes[strRep2][0], reverse, ref)
        #                 #     self.nodes[strRep2][0].addIn(node1, reverse, ref)
        #                 # print(strRep1, strRep2)
        #             else:
        #                 self.nodes[strRep1][0].addReadPos(pos)
        #                 self.nodes[strRep2][0].addReadPos(pos + 1)
        #         else:
        #             possibleCandidates = []
        #             for i in self.nodes[strRep1]:
        #                 for j in self.nodes[strRep2]:
        #                     if i in j.ins:
        #                         possibleCandidates.append([i,j])
        #             if len(possibleCandidates) == 1:
        #                 possibleCandidates[0][0].addOut(self.nodes[strRep2][0], reverse, ref)
        #                 possibleCandidates[0][1].addIn(self.nodes[strRep1][0],reverse, ref)
        #             else:
        #                 "uhhhadsfjlsd"
        #


    def topologicalOrder(self):
        self.cleanNodes()

        queue = []

        numberNodes = 0
        for node in self.nodes.values():
            for n in node:
                numberNodes += 1
                n.inDegree = len(n.ins)
                if len(n.ins) == 0:
                    queue.append(n)

        topOrder = []
        while queue:
            u = queue.pop()
            topOrder.append(u)

            for out in u.out:
                out.inDegree -= 1
                # If in-degree becomes zero, add it to queue
                if out.inDegree == 0:
                    queue.append(out)

        if len(topOrder) != numberNodes:
            return False
        else:
            return True






    def pruneGraph(self, minWeight):
        """
        Remove all edges with weight less than minWeight
        """
        allNodes = list(self.nodes.keys())
        for bunch in allNodes:
            for node in list(self.nodes[bunch]):
                #check all outs of this node
                outs = list(node.out.keys())
                for outNode in outs:
                    if node.ref and outNode.ref:
                        continue
                    elif sum(node.out[outNode]) < minWeight:
                        del node.out[outNode]
                        del outNode.ins[node]

        startingNodes = []
        for bunch in allNodes:
            for node in list(self.nodes[bunch]):
                if len(node.ins.keys()) == 0:
                    if len(node.out.keys()) == 0:
                        self.nodes[bunch].remove(node)
                    else:
                        startingNodes.append(node)



            if len(self.nodes[bunch]) == 0:
                del self.nodes[bunch]
        return startingNodes


    #     # remove all nodes that aren't connected by any edges
    #     self.removeNonConnectedNodes()
    #
    # def removeNonConnectedNodes(self):
    #     """
    #     Remove non-ref nodes that don't have ins or outs
    #     """
    #     numberNodes = len(self.nodes.keys())
    #     changed = True
    #     while changed:
    #         allNodes = list(self.nodes.keys())
    #         for node in allNodes:
    #             outs = list(self.nodes[node].out.keys())
    #             ins = list(self.nodes[node].ins.keys())
    #             if len(outs) <= 0 and self.nodes[node].ref is False:
    #                 self.removeNode(node)
    #             elif len(ins) <= 0 and self.nodes[node].ref is False:
    #                 self.removeNode(node)
    #         if numberNodes > len(self.nodes.keys()):
    #             numberNodes = len(self.nodes.keys())
    #         else:
    #             changed = False

    def removeNode(self, nodeName):
        """
        Remove a node from the graph with the given nodeName
        """
        node = self.nodes[nodeName]
        for ins in node.ins.keys():
            del self.nodes[ins].out[nodeName]
        for out in node.out.keys():
            del self.nodes[out].ins[nodeName]
        del self.nodes[nodeName]

    def cleanNodes(self):
        """
        dfs helper function, sets all nodes to unvisited state
        """
        for node in self.nodes.values():
            for n in node:
                n.visited = False

    # def search(self, node, discovered):
    #     discovered.append(node)
    #
    #     if not node.out.keys():  # no keys
    #         self.discovered.append(discovered)
    #
    #     for nextNode in node.out.keys():
    #         if nextNode not in discovered:
    #             self.search(nextNode, discovered.copy())
    #         else: # broke cycle, so add this path
    #             self.discovered.append(discovered)

    def search(self, node, discovered, pos):

        discovered[node].append(pos)

        if not node.out.keys():  # no keys
            if self.minMultiplicity(discovered):
                orderedDiscovered = [0]*(pos+1)
                for key,value in discovered.items():
                    for v in value:
                        orderedDiscovered[v] = key
                # print(mergeNodes(orderedDiscovered))
                self.discovered.append(orderedDiscovered)

        for nextNode in node.out.keys():
            if len(discovered[nextNode]) < nextNode.multiplicity[1]:
                copied = defaultdict(list)
                for key, valuelist in discovered.items():
                    copied[key] = valuelist.copy()
                self.search(nextNode, copied, pos + 1)
            else: # broke cycle, so add this path
                if self.minMultiplicity(discovered):
                    orderedDiscovered = [0] * (pos+1)
                    for key, value in discovered.items():
                        for v in value:
                            orderedDiscovered[v] = key
                    # print(mergeNodes(orderedDiscovered))
                    self.discovered.append(orderedDiscovered)


    def minMultiplicity(self, kmers):
        for key, value in self.nodes.items():
            for v in value:
                if len(kmers[v]) >= v.multiplicity[0] and len(kmers[v]) <= v.multiplicity[1]:
                    continue
                else:
                    return False
        return True

    def performSearch(self, startingNodes):
        self.discovered = []
        for node in startingNodes:
            self.search(node, defaultdict(list), 0)

    def findPaths(self, selectedPaths=False):
        chosenPaths = []
        if selectedPaths:
            for p in selectedPaths:
                chosenPaths.append(self.discovered[p])
        else:
            chosenPaths = self.discovered
        for path in chosenPaths:
            startRef = 0
            switch = False
            endRef = len(path)
            lastRef = 0
            switchRef = False
            for i in range(0, len(path)):
                # print(i, path[i].name)
                if path[i].ref:
                    if lastRef != 0 or path[i-1].ref is not True:
                        if lastRef + 1 in path[i].refPos:
                            lastRef += 1
                            if switchRef and endRef > i:
                                endRef = i + 1
                        else:
                            if switchRef == False:
                                switchRef = True
                                startRef = i
                            else:
                                lastRef += 1

                    else:
                        lastRef = path[i].refPos[0]

            # for i in range(0, len(path)):
            #     if path[i].ref:
            #         if switch and endRef > i: #works if there is a nonref node. TODO fix if ref to ref
            #             endRef = i + 1
            #         elif not switch and startRef < i:
            #             startRef = i
            #     else:
            #         switch = True
            if switch or switchRef:
                seq = mergeNodes(path[startRef:endRef])
                # should always start on ref
                ref1Pos = path[startRef].refPos[0]
                # what if don't end on ref? Ignore path.
                if path[endRef-1].ref:
                    ref2Pos = path[endRef-1].refPos[0]

                    #what if ref2 is before ref1? swap ref1 and ref2
                    swapped = False
                    if ref1Pos > ref2Pos:
                        ref1Pos, ref2Pos = ref2Pos, ref1Pos
                        swapped = True

                    ref = self.findRefFromPos(ref1Pos, ref2Pos)
                    minRef, minSeq, minPos = self.findMinRepresentation(ref, seq, ref1Pos, ref2Pos, swapped)
                    minCount, maxCount, average, fractionForward = self.nodeStats(path[startRef:endRef])
                    yield (self.chr, minPos, ".", minRef, minSeq, minCount, fractionForward)

    # def findPaths(self, selectedPaths=False):
    #     if selectedPaths == False:
    #         selectedPaths == [i for i in range(0, len(self.discovered))]
    #
    #     for i in range(0, len(self.discovered)):
    #         if i in selectedPaths:
    #             yield self.compare(self.discovered[i], self.discovered[self.refPath])
    #
    #
    # def compare(self, path1, path2):
    #     print('eh')
    #     if self.nodes[path1[0]].ref:
    #         minRef, minSeq, minPos = self.findMinRepresentation(path1, path2, self.nodes[path1[0]].refPos[0])
    #     elif self.nodes[path1[-1]].ref:
    #         minRef, minSeq, minPos = self.findMinRepresentation(path1, path2, self.nodes[path1[0]].refPos[-1], True)
    #     else:
    #         for p in range(0, len(path1)):
    #             if self.nodes[path1[p]].ref:
    #                 for r in range(0, len(path2)):
    #                     if path2[r] == path1[p]:
    #                         minRef, minSeq, minPos = self.findMinRepresentation(path1[:p], path2[:r], self.nodes[self.discovered[0]].refPos[-1], True)
    #                         if len(minSeq) == 0:
    #                             minRef, minSeq, minPos = self.findMinRepresentation(path1[p:], path2[r:], self.nodes[self.discovered[0]].refPos[-1])

        # minCount, maxCount, average, fractionForward = self.nodeStats(path2[startRef:endRef])
        # yield (self.chr, minPos, ".", minRef, minSeq, minCount, fractionForward)

        # yield (self.chr, minPos, ".", minRef, minSeq, 1, 1)




    def findMinRepresentation(self, ref, seq, ref1Pos, ref2Pos, swapped=False):
        """
        Reduce ref and seq sequences into the minimum representation seen in the vcf.
        """
        # if len(ref) == len(seq) and swapped is False: #SNP
        #     return ref[self.k:-self.k], seq[self.k:-self.k], ref1Pos+self.k
        # else:
        # thanks: https://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/
        while (min(len(seq), len(ref)) > 1 and seq[-1] == ref[-1]):
            seq = seq[:-1]
            ref = ref[:-1]
            # strip off identical prefixes and increment position
        while (min(len(seq), len(ref)) > 1 and seq[0] == ref[0]):
            seq = seq[1:]
            ref = ref[1:]
            ref1Pos += 1
        return ref, seq, ref1Pos

    def multiplicity(self, reads):
        for kmer in self.nodes.keys():
            countList = []
            for read in reads:
                countList.append(occurrences(read, kmer))

            self.nodes[kmer][0].multiplicity = [min(countList), max(countList)]
            if self.nodes[kmer][0].ref:
                self.nodes[kmer][0].refMultiplicity = occurrences(self.ref, kmer)



    def findCandidateStarts(self):
        """
        Find nodes with no ins.
        """
        candidateStarts = []
        for node in self.nodes.values():
            if len(node.ins) == 0:
                candidateStarts.append(node)
        return candidateStarts

    def nodeStats(self, kmers):
        """
        Find the minimum, maximum, and average counts of the given kmers
        """
        counts = []
        numberForward = 0
        for i in range(1, len(kmers)):
            counts.append(sum(kmers[i-1].out[kmers[i]]))
            numberForward += kmers[i-1].out[kmers[i]][0]
        if len(counts):
            minCount = min(counts)
            maxCount = max(counts)
            average =  sum(counts)/len(counts)


        return (minCount, maxCount, average, numberForward/sum(counts))

    def findRefFromPos(self, pos1, pos2):
        """
        Return the reference sequence given the two positions
        """
        start = pos1 - self.pos1
        end = pos2 - self.pos1 + self.k
        return self.ref[start:end]

    # def findVariants(self, best2Comb, bestPathForEachRead, best2AlignmentForEachRead, bestSumCount, variantDict, variants4eachPath):
    #     for c in range(0, len(best2Comb)):
    #         for CHROM, POS, ID, REF, ALT, MINCOUNT, FRACFORWARD in self.findPaths([best2Comb[c]]):
    #             myVariant = Variant(CHROM, POS, ID, REF, ALT, MINCOUNT,FRACFORWARD)
    #             if len(best2Comb) == 2:
    #                 FA = len(bestPathForEachRead[c]) / len(best2AlignmentForEachRead)
    #                 myVariant.FA = [FA]
    #                 AS = sum(best2AlignmentForEachRead) - bestSumCount
    #                 myVariant.AS = [AS]
    #             else:
    #                 FA = -1
    #                 myVariant.FA = [FA]
    #                 AS = -1
    #                 myVariant.AS = [AS]
    #             if POS in variantDict:
    #                 sameVariant = False
    #                 for variant in variantDict[POS]:
    #                     if variant.REF == REF and variant.ALT == ALT:
    #                         variant.NS = max(MINCOUNT, variant.NS)
    #                         variant.NW += 1
    #                         sameVariant = True
    #                         variant.FA.append(FA)
    #                         variant.AS.append(AS)
    #                         variant.FF.append(FRACFORWARD)
    #                 if sameVariant == False:
    #                     variantDict[POS].append(myVariant)
    #             else:
    #                 variantDict[POS] = [myVariant]
    #     return variantDict

    def findVariants(self, best2Comb, bestPathForEachRead, best2AlignmentForEachRead, bestSumCount, variantDict, variants4eachPath):
        for c in range(0, len(best2Comb)):
            for varList in variants4eachPath[best2Comb[c]].values():
                for myVariant in varList:
                    if len(best2Comb) == 2:
                        FA = len(bestPathForEachRead[c]) / len(
                            best2AlignmentForEachRead)
                        myVariant.FA = [FA]
                        AS = sum(best2AlignmentForEachRead) - bestSumCount
                        myVariant.AS = [AS]
                    else:
                        FA = -1
                        myVariant.FA = [FA]
                        AS = -1
                        myVariant.AS = [AS]
                    if myVariant.POS in variantDict:
                        sameVariant = False
                        for variant in variantDict[myVariant.POS]:
                            if variant.REF == myVariant.REF and variant.ALT == myVariant.ALT:
                                variant.NS = max(myVariant.NS, variant.NS)
                                variant.NW += 1
                                sameVariant = True
                                variant.FA.append(FA)
                                variant.AS.append(AS)
                                variant.FF.append(myVariant.FF[0])
                        if sameVariant == False:
                            variantDict[myVariant.POS].append(myVariant)
                    else:
                        variantDict[myVariant.POS] = [myVariant]
        return variantDict


    def storeVariants(self, POS, REF, ALT, path, variantDict):
        MINCOUNT, maxCount, average, FRACFORWARD = self.nodeStats(path)
        myVariant = Variant(self.chr, POS, ".", REF, ALT, MINCOUNT, FRACFORWARD)

        FA = -1
        myVariant.FA = [FA]
        AS = -1
        myVariant.AS = [AS]

        if POS in variantDict:
            sameVariant = False
            for variant in variantDict[POS]:
                if variant.REF == REF and variant.ALT == ALT:
                    variant.NS = max(MINCOUNT, variant.NS)
                    variant.NW += 1
                    sameVariant = True
                    # variant.FA.append(FA)
                    # variant.AS.append(AS)
                    variant.FF.append(FRACFORWARD)
            if sameVariant == False:
                variantDict[POS].append(myVariant)
        else:
            variantDict[POS] = [myVariant]
        return variantDict

    def printVariants(self, variantDict, index=False, debug=False):
        # print variants with pos smaller than start index
        for i in sorted(variantDict):
            if index is False or i < index:
                for variant in variantDict[i]:
                    array2print = [variant.CHROM, variant.POS, variant.ID,
                              variant.REF, variant.ALT, "50", "PASS",
                              "NS=" + str(variant.NS) + ";FF=" + str(
                                  sum(variant.FF) / len(variant.FF))
                              + ";FA=" + str(sum(variant.FA) / len(
                                  variant.FA)) + ";AS=" + str(
                                  sum(variant.AS) / len(
                                      variant.AS)) + ";NW=" + str(variant.NW)]

                    string2print = '\t'.join([str(i) for i in array2print] )

                    if debug:
                        logging.debug(string2print)
                    else:
                        print(string2print)
                del variantDict[i]
            else:
                break
        return variantDict

    def printGraph(self, cutoff=0):
        """
        Create an image representation of the graph.
        Can set a cutoff larger than pruning, but computations are run on pruned graph
        """
        dot = Digraph(comment=self.name, engine="dot")
        dot.attr(rankdir='LR') #draw from Left to Right
        for key, value in self.nodes.items():
            for n in value:
                if n.ref:
                    dot.node(str(n), key +"\n"+str(n.refPos) +"\n"+str(n.multiplicity), style='filled', fillcolor="#E1ECF4")
                else:
                    dot.node(str(n), key+"\n"+str(n.multiplicity))
                    for node, count in n.out.items():
                        if sum(count) > cutoff:
                            if node.ref == False:
                                dot.node(str(node), key+"\n"+str(n.multiplicity))
        for key, value in self.nodes.items():
            for n in value:
                for node, count in n.out.items():
                    if n.ref and node.ref:
                        dot.edge(str(n), str(node), xlabel=str(count))
                    elif sum(count) > cutoff:
                        dot.edge(str(n), str(node), xlabel=str(count))
        # print(dot.source)
        dot.render('test-output/' + self.name, view=True)


# def deBruijn(graph, k, text, reverse=True, ref=False, final=False):
#     """
#     Construct a deBruijn graph from a string.
#     """
#     nodes = [text[0:0 + k]]
#     for i in range(0, len(text) - k):
#         node1 = text[i:i + k]
#         node2 = text[i + 1:i + k + 1]
#         if ref and final is False:
#             if node1 in graph.nodes and node2 in graph.nodes:
#                 return False
#         graph.addEdge(node1, node2, reverse, ref)
#
#         nodes.append(node2)
#     return nodes

def deBruijn(graph, k, text, firstPos=0, reverse=True, ref=False):
    """
    Construct a deBruijn graph from a string.
    """
    pos = 0
    myNode = Node(text[0:0 + k], ref)
    if ref:
        myNode.addRefPos(pos + firstPos)
    nodes = [myNode]
    for i in range(1, len(text) - k+1):
        strRep = text[i:i + k]
        myNode = Node(strRep, ref)
        if ref:
            pos += 1
            myNode.addRefPos(pos + firstPos)
        graph.addEdge(nodes[-1].name, strRep,nodes[-1], myNode, reverse, ref)
        nodes.append(myNode)
    return nodes


def occurrences(string, sub):
    # https://stackoverflow.com/questions/2970520/string-count-with-overlapping-occurrences
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count

def mergeNodes(kmers):
    """
    Create the string formed by kmers. Kmers in correct order.
    """
    text = kmers[0].name
    for i in range(1, len(kmers)):
        text = text + kmers[i].name[-1]
    return text

def findDifference(ref, seq, ref1Pos, path):
    # print(ref,"\n"+seq, ref1Pos, path)
    savedRef = ref
    savedSeq = seq
    savedPath = path.copy()
    savedPathPointer = 0
    pointer = 0
    while (min(len(seq), len(ref)) > 1):
        if seq[0] == ref[0]:
            seq = seq[1:]
            ref = ref[1:]
            ref1Pos += 1
            pointer += 1
            path = path[1:]
            savedPathPointer += 1

        else:
            increment = 0
            while seq[increment] == "-" or ref[increment] == "-" or ref[increment] != seq[increment] :
                increment += 1
                if pointer + increment >= len(savedSeq):
                    break
            if increment:
                print(ref, seq)
                if ref[0] == "-":
                    path2return = howFar(savedPath, savedPathPointer-1, savedPathPointer+increment)
                    if path2return is not False:
                        # print("if", savedRef[pointer-1:pointer+increment].replace("-", ""), savedSeq[pointer-1:pointer+increment].replace("-", ""))
                        yield(savedRef[pointer-1:pointer+increment].replace("-", ""), savedSeq[pointer-1:pointer+increment].replace("-", ""), ref1Pos-1, path2return)
                    savedPathPointer += increment
                    path = path[increment:]
                elif seq[0] == "-":
                    path2return = howFar(savedPath, savedPathPointer, savedPathPointer+increment+1)
                    if path2return is not False:
                        # print("elif",savedRef[pointer:pointer + increment+1].replace("-", ""), savedSeq[pointer:pointer + increment+1].replace("-", "") )
                        yield(savedRef[pointer:pointer + increment+1].replace("-", ""), savedSeq[pointer:pointer + increment+1].replace("-", ""), ref1Pos-1, path2return)
                else:
                    path2return = howFar(savedPath, savedPathPointer,  savedPathPointer + increment)
                    if path2return is not False:
                        # print("else", ref[:increment].replace("-", ""), seq[:increment].replace("-", ""))
                        yield(ref[:increment].replace("-", ""), seq[:increment].replace("-", ""), ref1Pos,path2return)
                    savedPathPointer += increment
                    path = path[increment:]
                seq = seq[increment:]
                ref = ref[increment:]
                pointer += increment
                ref1Pos += increment
            else:
                seq = seq[1:]
                ref = ref[1:]
                path = path[1:]
                ref1Pos += 1
                pointer += 1
                savedPathPointer += 1

def howFar(savedPath, begin, end, k=8):
    if end >= len(savedPath):
        return False
    if begin > end:
        return False

    if savedPath[begin].ref and savedPath[end].ref:
        return savedPath[begin: end]

    newBegin = begin
    newEnd = end
    # if savedPath[begin].ref:
    while newEnd < len(savedPath) and savedPath[newEnd].ref == False and newEnd - end < k-1:
        newEnd += 1
    # if savedPath[end].ref:
    while savedPath[newBegin].ref == False and  begin - newBegin < k-1 and newBegin > 0:
        newBegin -= 1

    return savedPath[newBegin: newEnd]



def argParser():
    """
    Parse user inputs.
    """
    parser=argparse.ArgumentParser(add_help=True)
    parser.add_argument("--outputFile", "-o",
                        type=str,
                        help="specify the file name of the output file")
    parser.add_argument("--inputBam", "-i",
                        type=str,
                        help="specify the file name of the input BAM file")
    parser.add_argument("--inputRef", "-r",
                        type=str,
                        help="specify the file name of the input reference file")

    return vars(parser.parse_args())

def main():
    """
    Collect user input and find variants.
    """
    args = argParser()
    parseBam(args["inputBam"], args["inputRef"], args["outputFile"])



if __name__ == "__main__":
    main()