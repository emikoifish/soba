# python3 parser.py -i HG002_PacBio_GRCh38.bam -o test.txt

import pysam
import argparse
import sys
from graphviz import Digraph
from pyfasta import Fasta
import math

def parseBam(inputFile, outputFile):
    samfile = pysam.AlignmentFile(inputFile, "rb")


    chr = "chr20"
    startIndex = 1120000
    endIndex = 1130000
    # startIndex = 1000000
    # endIndex = 2000000
    k = 8
    windowSize = 30
    windowOverlap = 5

    currentPositionVariants = {}

    if endIndex - startIndex < windowSize:
        endIndex = startIndex + windowSize
        print("endIndex needs to be larger than or equal to startIndex + windowSize. Setting endIndex to: ", endIndex)

    print("CHROM\tPOS\tID\tREF\tALT\tMINCOUNT")
    for interval in range(startIndex, endIndex - windowSize+1, windowOverlap):
        index1 = interval
        index2 = interval + windowSize

        myGraph = Graph(chr+"_"+str(index1)+"_"+str(index2), chr, index1, index2, k)
        index1 -= 1

        # get the ref sequence
        f = Fasta('GCA_000001405.15_GRCh38_no_alt_analysis_set.fna')
        ref = f["chr20  AC:CM000682.2  gi:568336004  LN:64444167  rl:Chromosome  M5:b18e6c531b0bd70e949a7fc20859cb01  AS:GRCh38"][index1:index2]

        refNodes = deBruijn(myGraph, k, ref, ref=True)
        i = index1 + 1
        for node in refNodes:
            myGraph.nodes[node].ref = True
            myGraph.nodes[node].addRefPos(i)
            i += 1
        myGraph.ref = ref

        for read in samfile.fetch('chr20', index1, index2):
            refPos = read.get_reference_positions(full_length=True)
            # cigar = read.cigartuples
            # print(cigar)

            if index1 in refPos:
                seqFirstIndex = refPos.index(index1)
            else:
                seqFirstIndex = findClosestIndex(refPos, index1, True)

            if index2 in refPos:
                seqLastIndex = refPos.index(index2)
            else:
                seqLastIndex = findClosestIndex(refPos, index2, False)

            windowedRead = read.query_alignment_sequence[seqFirstIndex:seqLastIndex]

            deBruijn(myGraph, k, windowedRead, read.is_reverse)

        myGraph.pruneGraph(5)
        # myGraph.collapsedGraph()
        
        # find the variants and add to currentPositionVariants
        for CHROM, POS, ID, REF, ALT, MINCOUNT in myGraph.findVariants(refNodes):
            if POS in currentPositionVariants:
                sameVariant = False
                for variant in currentPositionVariants[POS]:
                    if variant[2] == REF and variant[3] == ALT:
                        variant[4] = max(MINCOUNT, variant[4])
                        sameVariant = True
                if sameVariant == False:
                    currentPositionVariants[POS].append([CHROM, ID, REF, ALT, MINCOUNT])
            else:
                currentPositionVariants[POS] = [[CHROM, ID, REF, ALT, MINCOUNT]]

        # print variants with pos smaller than start index
        for i in sorted(currentPositionVariants):
            if i < index1:
                for variant in currentPositionVariants[i]:
                    print(variant[0], i, variant[1],variant[2],variant[3],variant[4], sep="\t")
                del currentPositionVariants[i]
            else:
                break

        # myGraph.printGraph()

    # print the rest of the variants that haven't been removed from the list
    for i in sorted(currentPositionVariants):
        for variant in currentPositionVariants[i]:
            print(variant[0], i, variant[1], variant[2], variant[3], variant[4], sep="\t")
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
                    if i == 0 or lastRealRefPos==0:
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
            self.ins[newIn] = [0,0]
            self.ins[newIn][reverse] += addCount

    def addOut(self, newOut, reverse, ref):
        """Add a new out and in to node, if haven't seen out before, add to dictionary."""
        addCount = 1
        if ref:
            addCount = 0
        if newOut in self.out:
            self.out[newOut][reverse] += addCount
        else:
            self.out[newOut] = [0,0]
            self.out[newOut][reverse] += addCount

    def addRefPos(self, pos):
        self.refPos.append(pos)

class Graph:
    """
    A collection of nodes.
    """
    def __init__(self, name,chr, pos1, pos2, k):
        self.name = name
        self.nodes = {}
        self.totalEdges = 0
        self.ref = ""
        self.chr = chr
        self.pos1 = pos1
        self.pos2 = pos2
        self.k = k

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

        # remove all nodes that aren't connected by any edges
        self.removeNonConnectedNodes()

    def removeNonConnectedNodes(self):
        """
        Remove non-ref nodes that don't have ins or outs
        """
        numberNodes = len(self.nodes.keys())
        changed = True
        while changed:
            allNodes = list(self.nodes.keys())
            for node in allNodes:
                outs = list(self.nodes[node].out.keys())
                ins = list(self.nodes[node].ins.keys())
                if len(outs) <= 0 and self.nodes[node].ref == False:
                    self.removeNode(node)
                elif len(ins) <= 0 and self.nodes[node].ref == False:
                    self.removeNode(node)
            if numberNodes > len(self.nodes.keys()):
                numberNodes = len(self.nodes.keys())
            else:
                changed = False

    def removeNode(self, nodeName):
        node = self.nodes[nodeName]
        for ins in node.ins.keys():
            del self.nodes[ins].out[nodeName]
        for out in node.out.keys():
            del self.nodes[out].ins[nodeName]
        del self.nodes[nodeName]

    def dfs(self, start):
        """
        Do a Depth First Search based search of the graph for paths
        """
        finishedPaths = []
        stack = [(start, [])]
        while stack:
            node, prevList = stack.pop()
            self.nodes[node].visited = True
            prevList.append(node)
            for nextNode, count in self.nodes[node].out.items():
                if self.nodes[nextNode].ref == True:
                    finishedList = prevList.copy()
                    finishedList.append(nextNode)
                    finishedPaths.append(finishedList)
                elif self.nodes[nextNode].visited == False and self.nodes[nextNode].ref == False:
                    stack.append((nextNode, prevList.copy()))
        return finishedPaths

    def findVariants(self, refNodes):
        """
        Locate variants in the given graph.
        """
        # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
        for refNode in refNodes:
            for path in self.dfs(refNode):
                ref1 = path[0]
                ref2 = path[-1]

                ref1Pos = self.nodes[ref1].refPos
                ref2Pos = self.nodes[ref2].refPos
                if len(path) > 2:
                    seq = self.mergeNodes(path)
                    if ref1Pos[0] < ref2Pos[0]:
                        ref = self.findRefFromPos(ref1Pos[0],ref2Pos[0])

                        minRef, minSeq, minPos = self.findMinRepresentation(ref, seq, ref1Pos[0],ref2Pos[0])

                        yield (self.chr, minPos, "ID",minRef, minSeq, self.nodeStats(path)[0])
                        # print(ref1Pos, ref2Pos, ref, seq)
                    else:
                        # cycle in graph, do not print
                        #print("Cycle in graph")
                        continue

                # else:
                #     # if deletions this should be called
                #     # however, make sure that counts are high to be called?
                #     # actually more thought should be put into this
                #     nextToEachOther = None
                #     for pos1 in ref1Pos:
                #         for pos2 in ref2Pos:
                #             if nextToEachOther:
                #                 break
                #             elif abs(pos2 - pos1) > 1:
                #                 nextToEachOther = False
                #                 p1 = pos1
                #                 p2 = pos2
                #             else:
                #                 nextToEachOther = True
                #                 break
                #     if nextToEachOther == False:
                #         if self.nodes[ref1].out[ref2] > 1:
                #             print(path, ref1Pos, ref2Pos, ref1+ref2[-1])



    def findMinRepresentation(self, ref, seq, ref1Pos, ref2Pos):
        """
        Reduce ref and seq sequences into the minimum representation seen in the vcf.
        """
        if len(ref) == len(seq): #SNP
            return ref[self.k:-self.k], seq[self.k:-self.k], ref1Pos+self.k
        else:
            # thanks: https://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/
            while ( min(len(seq), len(ref)) > 1 and seq[-1] == ref[-1]):
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

    def mergeNodes(self, kmers):
        """
        Create the string formed by kmers. Kmers in correct order.
        """
        text = kmers[0]
        for i in range(1, len(kmers)):
            text = text + kmers[i][-1]
        return text

    def nodeStats(self, kmers):
        """
        Find the minimum, maximum, and average counts of the given kmers
        """
        counts = []
        for i in range(1, len(kmers)):
            counts.append(sum(self.nodes[kmers[i-1]].out[kmers[i]]))
        minCount = min(counts)
        maxCount = max(counts)
        average = sum(counts)/len(counts)
        return (minCount, maxCount, average)

    def findRefFromPos(self, pos1, pos2):
        """
        Return the reference sequence given the two positions
        """
        start = pos1 - self.pos1
        end = pos2 - self.pos1 + self.k
        return self.ref[start:end]

    # Logic to merge nodes into larger nodes
    # def replaceNodes(self, mergeNodes, replacementNode, mergedCounts, ref=False):
    #     self.nodes[replacementNode] = Node(replacementNode, ref=ref, merged=True, mergedCounts=mergedCounts)
    #     self.nodes[replacementNode].ins = self.nodes[mergeNodes[0]].ins
    #     self.nodes[replacementNode].out = self.nodes[mergeNodes[-1]].out
    #
    #     for inNode, count in self.nodes[mergeNodes[0]].ins.items():
    #         del self.nodes[inNode].out[mergeNodes[0]]
    #         self.nodes[inNode].out[replacementNode] = min(mergedCounts)
    #
    #     for outNode, count in self.nodes[mergeNodes[-1]].out.items():
    #         del self.nodes[outNode].ins[mergeNodes[-1]]
    #         self.nodes[outNode].ins[replacementNode] = min(mergedCounts)
    #
    #     for node in mergeNodes:
    #         del self.nodes[node]
    #
    # def mergeable(self, currentNode):
    #     print(currentNode.name, currentNode.out.values(), currentNode.ins.values())
    #     if len(currentNode.out) == 1 and len(currentNode.ins) <= 1:
    #         if currentNode.ref == self.nodes[list(currentNode.out.keys())[0]].ref:
    #             for key, value in currentNode.out.items():
    #                 return value
    #     return False
    #
    # def collapsedGraph(self):
    #     tempNodes = self.findCandidateStarts()
    #     while tempNodes: # for each tree, should only be 1
    #         currentNode = tempNodes.pop()
    #         mergeableNodes = []
    #         toSearch = [currentNode]
    #         mergedCounts = []
    #         while len(toSearch): # traverse down tree
    #             currentNode.merged = True
    #             count = self.mergeable(currentNode)
    #             if count:
    #                 mergeableNodes.append(currentNode)
    #                 mergedCounts.append(count)
    #                 currentNode=self.nodes[list(currentNode.out.keys())[0]]
    #
    #             else:
    #                 if len(currentNode.out) > 1:
    #                     mergeableNodes.append(currentNode)
    #                     mergedCounts.append(count)
    #                 if len(mergeableNodes) > 0:
    #                     print("merging")
    #
    #                     mergedNodes = self.mergeNodes([i.name for i in mergeableNodes])
    #                     print(mergedNodes)
    #                     print(mergedCounts)
    #                     self.replaceNodes([i.name for i in mergeableNodes], mergedNodes, mergedCounts)
    #                     mergeableNodes = []
    #                     mergedCounts = []
    #                 for key in currentNode.out.keys():
    #                     if self.nodes[key].merged == False:
    #                         toSearch.append(self.nodes[key])
    #                 currentNode = toSearch.pop()


    def printGraph(self, cutoff=0):
        """
        Create an image representation of the graph.
        Can set a cutoff larger than pruning, but computations are run on pruned graph
        """
        dot = Digraph(comment=self.name, engine="dot")
        for key, value in self.nodes.items():
            if value.ref:
                dot.node(key, key, style='filled', fillcolor="#E1ECF4")
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


def deBruijn(graph, k, text, reverse=True, ref=False):
    """
    Construct a deBruijn graph from a string.
    """
    nodes = [text[0:0 + k]]
    for i in range(0, len(text) - k):
        node1 = text[i:i + k]
        node2 = text[i + 1:i + k + 1]
        graph.addEdge(node1, node2, reverse, ref)
        nodes.append(node2)
    return nodes


def argParser():
    parser=argparse.ArgumentParser(add_help=True)
    parser.add_argument("--outputFile","-o",
                        type=str,
                        help="specify the file name of the output file")
    parser.add_argument("--inputFile","-i",
                        type=str,
                        help="specify the file name of the input file")
    return vars(parser.parse_args())

def main():
    args = argParser()
    parseBam(args["inputFile"],args["outputFile"])



if __name__ == "__main__":
    main()