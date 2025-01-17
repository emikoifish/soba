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
    startIndex = 2145400
    endIndex = startIndex + 50
    # startIndex = 2145400
    # endIndex = startIndex + 10000
    # startIndex = 2000000
    # endIndex = 2023000

    # startIndex = 2000000
    # endIndex = 3000000
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
                refNodes = deBruijn(myGraph, changedk, ref, ref=True)
                while refNodes == False and changedk <= 15:
                    changedk += 1
                    myGraph = Graph(chr + "_" + str(index1) + "_" + str(index2), chr, index1, index2, changedk)
                    refNodes = deBruijn(myGraph, changedk, ref, ref=True)
                if refNodes == False:
                    refNodes = deBruijn(myGraph, changedk, ref, ref=True, final=True)

                for node in refNodes:
                    myGraph.nodes[node].ref = True
                    myGraph.nodes[node].addRefPos(i)
                    i += 1
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
                        if len(windowedRead):
                            savedWindowedReads.append(windowedRead)
                        # print(refPosStart, refPosEnd, windowedRead)
                        deBruijn(myGraph, changedk, windowedRead, read.is_reverse)
                    else:
                        # remove from memory
                        removeIndexes.append(i)
                for index in removeIndexes[::-1]:
                    del cursors[index]


                myGraph.pruneGraph(2)

                myGraph.printGraph()

                # find paths
                myGraph.performSearch(refNodes[0])

                # make sure to add the ref seq if doesn't exist.
                refPath = -1
                for i in range(0, len(myGraph.discovered)):
                    if myGraph.discovered[i] == refNodes:
                        refPath = i
                if refPath == -1:
                    myGraph.discovered.append(refNodes)
                    refPath = len(myGraph.discovered) - 1



                alignments = []
                # which path do you belong to the best?
                for path in myGraph.discovered:
                    newAlignments = []
                    for read in savedWindowedReads:
                        newAlignments.append(max(i[2] for i in pairwise2.align.globalms(mergeNodes(path), read, 2, -1, -.5, -.1)))
                    alignments.append(newAlignments)

                # instead of finding the best path, we want to find the reference path.
                refSumCount = sum(alignments[refPath])
                outFile.write("position: "+ str(index1) + "\t numPaths" + str(len(myGraph.discovered))+"\n")
                if len(myGraph.discovered) > 1:
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
                        outFile.write("3\n")


                    # find the variants and add to currentPositionVariants
                    variantDict = myGraph.findVariants(best2Comb, bestPathForEachRead, best2AlignmentForEachRead, refSumCount, variantDict)
                    variantDict = myGraph.printVariants(variantDict, index1)

                    allVariantsDict = myGraph.findVariants([i for i in range(0, len(myGraph.discovered))],bestPathForEachRead, best2AlignmentForEachRead,refSumCount, allVariantsDict)
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
    def __init__(self, name, ref=False, merged=False, mergedCounts=[]):
        """Store the node attributes and add node to dictionary of nodes."""
        self.out = {}
        self.ins = {}
        self.name = name
        self.ref = ref
        self.visited = False
        self.refPos = []
        self.merged = merged
        self.mergedCounts = mergedCounts

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

    def addRefPos(self, pos):
        """
        Add a reference pos to a ref node
        """
        self.refPos.append(pos)

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
        self.nodes = {}
        self.totalEdges = 0
        self.ref = ""
        self.chr = chr
        self.pos1 = pos1
        self.pos2 = pos2
        self.k = k
        self.discovered = []

    def addEdge(self, node1, node2, reverse, ref):
        """
        Add an edge from node1 to node2
        """
        if node1 not in self.nodes:
            self.nodes[node1] = Node(node1)
        if node2 not in self.nodes:
            self.nodes[node2] = Node(node2)
        self.nodes[node1].addOut(node2, reverse, ref)
        self.nodes[node2].addIn(node1, reverse, ref)

    def pruneGraph(self, minWeight):
        """
        Remove all edges with weight less than minWeight
        """
        allNodes = list(self.nodes.keys())
        for node in allNodes:

            #check all outs of this node
            outs = list(self.nodes[node].out.keys())
            for outNode in outs:
                # don't remove any ref nodes
                if self.nodes[node].ref and self.nodes[outNode].ref:
                    continue
                elif sum(self.nodes[node].out[outNode]) < minWeight:
                    del self.nodes[node].out[outNode]
                    del self.nodes[outNode].ins[node]

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
            node.visited = False

    def search(self, node, discovered):
        discovered.append(node)

        if not self.nodes[node].out.keys():  # no keys
            self.discovered.append(discovered)

        for nextNode in self.nodes[node].out.keys():
            if nextNode not in discovered:
                self.search(nextNode, discovered.copy())
            else: # broke cycle, so add this path
                self.discovered.append(discovered)

    def performSearch(self, startingNode):
        self.discovered = []
        self.search(startingNode, [])

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
            for i in range(0, len(path)):
                if self.nodes[path[i]].ref:
                    if switch and endRef > i: #works if there is a nonref node. TODO fix if ref to ref
                        endRef = i + 1
                    elif not switch and startRef < i:
                        startRef = i
                else:
                    switch = True
            if switch:
                seq = mergeNodes(path[startRef:endRef])
                # should always start on ref
                ref1Pos = self.nodes[path[startRef]].refPos[0]
                # what if don't end on ref? Ignore path.
                if self.nodes[path[endRef-1]].ref:
                    ref2Pos = self.nodes[path[endRef-1]].refPos[0]

                    #what if ref2 is before ref1? swap ref1 and ref2
                    swapped = False
                    if ref1Pos > ref2Pos:
                        ref1Pos, ref2Pos = ref2Pos, ref1Pos
                        swapped = True

                    ref = self.findRefFromPos(ref1Pos, ref2Pos)
                    minRef, minSeq, minPos = self.findMinRepresentation(ref, seq, ref1Pos, ref2Pos, swapped)
                    minCount, maxCount, average, fractionForward = self.nodeStats(path[startRef:endRef])
                    yield (self.chr, minPos, ".", minRef, minSeq, minCount, fractionForward)

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
            counts.append(sum(self.nodes[kmers[i-1]].out[kmers[i]]))
            numberForward += self.nodes[kmers[i-1]].out[kmers[i]][0]
        minCount = min(counts)
        maxCount = max(counts)
        average = sum(counts)/len(counts)

        return (minCount, maxCount, average, numberForward/sum(counts))

    def findRefFromPos(self, pos1, pos2):
        """
        Return the reference sequence given the two positions
        """
        start = pos1 - self.pos1
        end = pos2 - self.pos1 + self.k
        return self.ref[start:end]

    def findVariants(self, best2Comb, bestPathForEachRead, best2AlignmentForEachRead, bestSumCount, variantDict):
        for c in range(0, len(best2Comb)):
            for CHROM, POS, ID, REF, ALT, MINCOUNT, FRACFORWARD in self.findPaths([best2Comb[c]]):
                myVariant = Variant(CHROM, POS, ID, REF, ALT, MINCOUNT,FRACFORWARD)
                if len(best2Comb) == 2:
                    FA = len(bestPathForEachRead[c]) / len(best2AlignmentForEachRead)
                    myVariant.FA = [FA]
                    AS = sum(best2AlignmentForEachRead) - bestSumCount
                    myVariant.AS = [AS]
                else:
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
                            variant.FA.append(FA)
                            variant.AS.append(AS)
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
        for key, value in self.nodes.items():
            if value.ref:
                dot.node(key, key +"\n"+str(value.refPos), style='filled', fillcolor="#E1ECF4")
            else:
                for node, count in value.out.items():
                    if sum(count) > cutoff:
                        dot.node(key, key)
        for key, value in self.nodes.items():
            for node, count in value.out.items():
                if self.nodes[key].ref and self.nodes[node].ref:
                    dot.edge(key, node, xlabel=str(count))
                elif sum(count) > cutoff:
                    dot.edge(key, node, xlabel=str(count))
        # print(dot.source)
        dot.render('test-output/' + self.name, view=True)


def deBruijn(graph, k, text, reverse=True, ref=False, final=False):
    """
    Construct a deBruijn graph from a string.
    """
    nodes = [text[0:0 + k]]
    for i in range(0, len(text) - k):
        node1 = text[i:i + k]
        node2 = text[i + 1:i + k + 1]
        if ref and final is False:
            if node1 in graph.nodes and node2 in graph.nodes:
                return False
        graph.addEdge(node1, node2, reverse, ref)
        nodes.append(node2)
    return nodes

def mergeNodes(kmers):
    """
    Create the string formed by kmers. Kmers in correct order.
    """
    text = kmers[0]
    for i in range(1, len(kmers)):
        text = text + kmers[i][-1]
    return text

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