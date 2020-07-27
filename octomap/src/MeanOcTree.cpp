/*
 * OctoMap - An Efficient Probabilistic 3D Mapping Framework Based on Octrees
 * http://octomap.github.com/
 *
 * Copyright (c) 2009-2013, K.M. Wurm and A. Hornung, University of Freiburg
 * All rights reserved.
 * License: New BSD
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Freiburg nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <octomap/MeanOcTree.h>

namespace octomap {


  // node implementation  --------------------------------------
  std::ostream& MeanOcTreeNode::writeData(std::ostream &s) const {
    s.write((const char*) &value, sizeof(value)); // occupancy
    s.write((const char*) &mean, sizeof(double)); // mean
    s.write((const char*) &n_observations, sizeof(long)); // number of observations

    return s;
  }

  std::istream& MeanOcTreeNode::readData(std::istream &s) {
    s.read((char*) &value, sizeof(value)); // occupancy
    s.read((char*) &mean, sizeof(double)); // mean
    s.read((char*) &n_observations, sizeof(long)); // number of observations

    return s;
  }

  void MeanOcTreeNode::addObservation(double i)
  {
    mean = ((mean * n_observations) + i);
    mean /= ++n_observations;
  }

  double MeanOcTreeNode::getAverageChildMean() const {
    double sum_means = 0;
    long sum_observations = 0;

    if (children != NULL){
      for (int i=0; i<8; i++) {
        MeanOcTreeNode* child = static_cast<MeanOcTreeNode*>(children[i]);

        if (child != NULL && child->isObserved()) {
          sum_means += child->getMean();
          sum_observations += child->getObservations();
        }
      }
    }

    if (sum_observations > 0) {
      return sum_means / double(sum_observations);
    }
    else { // no child had a color other than white
      return 0.0;
    }
  }


  void MeanOcTreeNode::updateMeanChildren() {
    mean = getAverageChildMean();
  }


  // tree implementation  --------------------------------------
  MeanOcTree::MeanOcTree(double in_resolution)
  : OccupancyOcTreeBase<MeanOcTreeNode>(in_resolution) 
  {
    meanOcTreeMemberInit.ensureLinking();
    intensity_threshold_ = 0.0;
  };

  MeanOcTree::MeanOcTree(double in_resolution, double i_thresh)
  : OccupancyOcTreeBase<MeanOcTreeNode>(in_resolution)
  {
    meanOcTreeMemberInit.ensureLinking();
    intensity_threshold_ = i_thresh;
  }

  MeanOcTreeNode* MeanOcTree::setNodeMean(const OcTreeKey& key, double m, long n)
  {
    MeanOcTreeNode* node = search (key);
    if (node != 0) {
      node->setMean(m,n);
    }
    return node;
  }

  MeanOcTreeNode* MeanOcTree::setNodeMean(const OcTreeKey& key,
                                             double m) {
    MeanOcTreeNode* node = search (key);
    if (node != 0) {
      node->setMean(m);
    }
    return node;
  }

  MeanOcTreeNode* MeanOcTree::setNodeObservations(const OcTreeKey& key, long n)
  {
    MeanOcTreeNode* node = search (key);
    if (node != 0)
    {
      node->setObservations(n);
    }
    return node;
  }

  MeanOcTreeNode* MeanOcTree::addObservation(const OcTreeKey& key, double i)
  {
    MeanOcTreeNode* node = search (key);
    if (node != 0)
    {
      node->addObservation(i);
    }
    return node;
  }

  bool MeanOcTree::pruneNode(MeanOcTreeNode* node) {
    if (!isNodeCollapsible(node))
      return false;

    // set value to children's values (all assumed equal)
    node->copyData(*(getNodeChild(node, 0)));

    if (node->isObserved()) // TODO: Unsure how to handle n_observations here
      node->setMean(node->getAverageChildMean());

    // delete children
    for (unsigned int i=0;i<8;i++) {
      deleteNodeChild(node, i);
    }
    delete[] node->children;
    node->children = NULL;

    return true;
  }

  bool MeanOcTree::isNodeCollapsible(const MeanOcTreeNode* node) const{
    // all children must exist, must not have children of
    // their own and have the same occupancy probability
    if (!nodeChildExists(node, 0))
      return false;

    const MeanOcTreeNode* firstChild = getNodeChild(node, 0);
    if (nodeHasChildren(firstChild))
      return false;

    for (unsigned int i = 1; i<8; i++) 
    {
      // compare nodes only using their occupancy, ignoring color for pruning
      if (!nodeChildExists(node, i) 
          || nodeHasChildren(getNodeChild(node, i)) 
          || !(getNodeChild(node, i)->getValue() == firstChild->getValue()))
        return false;
    }

    return true;
  }



  void MeanOcTree::updateInnerOccupancy() {
    this->updateInnerOccupancyRecurs(this->root, 0);
  }

  void MeanOcTree::updateInnerOccupancyRecurs(MeanOcTreeNode* node, unsigned int depth) {
    // only recurse and update for inner nodes:
    if (nodeHasChildren(node)){
      // return early for last level:
      if (depth < this->tree_depth){
        for (unsigned int i=0; i<8; i++) {
          if (nodeChildExists(node, i)) {
            updateInnerOccupancyRecurs(getNodeChild(node, i), depth+1);
          }
        }
      }
      node->updateOccupancyChildren();
      node->updateMeanChildren();
    }
  }


  std::ostream& operator<<(std::ostream& out, MeanOcTreeNode const* node) {
    return out << std::fixed << std::setprecision(3) 
      << '(' << node->getMean() << ' ' << node->getObservations() << ')';
  }


  MeanOcTree::StaticMemberInitializer MeanOcTree::meanOcTreeMemberInit;

} // end namespace
