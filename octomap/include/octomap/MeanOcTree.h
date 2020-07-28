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

#ifndef OCTOMAP_MEAN_OCTREE_H
#define OCTOMAP_MEAN_OCTREE_H


#include <iostream>
#include <octomap/OcTreeNode.h>
#include <octomap/OccupancyOcTreeBase.h>
#include <iomanip>

namespace octomap {
  
  // forward declaraton for "friend"
  class MeanOcTree;
  
  // node definition
  class MeanOcTreeNode : public OcTreeNode {    
  public:
    friend class MeanOcTree; // needs access to node children (inherited)

    MeanOcTreeNode() : OcTreeNode() {}

    MeanOcTreeNode(const MeanOcTreeNode& rhs) 
      : OcTreeNode(rhs), mean(rhs.mean), n_observations(rhs.n_observations) {}

    bool operator==(const MeanOcTreeNode& rhs) const{
      return (rhs.mean == mean && rhs.n_observations == n_observations);
    }
    
    void copyData(const MeanOcTreeNode& from){
      OcTreeNode::copyData(from);
      this->mean =  from.getMean();
      this->n_observations = from.getObservations();
    }
        
    inline double getMean() const { return mean; }
    inline long getObservations() const {return n_observations; }
    inline void  setMean(double m) {this->mean = m; }
    inline void setMean(double m, long n) {this->mean = m; this->n_observations = n; }
    inline void setObservations(long n){this->n_observations = n; }


    // has the node been observed?
    inline bool isObserved() const { 
      return (n_observations != 0); 
    }

    void updateMeanChildren();

    void addObservation(double i); // adds observation of cell with intensity i
    void updateLogOdds(); // updates log-odds occupancy based on mean intensity

    double getAverageChildMean() const;
  
    // file I/O
    std::istream& readData(std::istream &s);
    std::ostream& writeData(std::ostream &s) const;
    
  protected:
    double mean;
    long n_observations;
  };


  // tree definition
  class MeanOcTree : public OccupancyOcTreeBase <MeanOcTreeNode> {

  public:
    /// Default constructor, sets resolution of leafs
    MeanOcTree(double resolution);

    MeanOcTree(double resolution, double i_thresh);

    void setIntensityThreshold(double i_thresh) { intensity_threshold_ = i_thresh; } 

    /// virtual constructor: creates a new object of same type
    /// (Covariant return type requires an up-to-date compiler)
    MeanOcTree* create() const {return new MeanOcTree(resolution); }

    std::string getTreeType() const {return "MeanOcTree";}
    
     /**
     * Prunes a node when it is collapsible. This overloaded
     * version only considers the node occupancy for pruning,
     * different colors of child nodes are ignored.
     * @return true if pruning was successful
     */
    virtual bool pruneNode(MeanOcTreeNode* node);
    
    virtual bool isNodeCollapsible(const MeanOcTreeNode* node) const;
       
    // set node mean and/or number of observations at given key or coordinate. 
    MeanOcTreeNode* setNodeMean(const OcTreeKey& key, double m);
    MeanOcTreeNode* setNodeMean(const OcTreeKey& key, double m, long n);
    MeanOcTreeNode* setNodeObservations(const OcTreeKey& key, long n);
    MeanOcTreeNode* addObservation(const OcTreeKey& key, double i);


    MeanOcTreeNode* setNodeMean(float x, float y, 
                                 float z, double m) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return setNodeMean(key,m);
    }

    MeanOcTreeNode* setNodeMean(float x, float y, float z, double m, long n)
    {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return setNodeMean(key,m,n);
    }

    MeanOcTreeNode* setNodeObservations(float x, float y, float z, long n)
    {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z),key)) return NULL;
      return setNodeObservations(key,n);
    }

    MeanOcTreeNode* addObservation(float x, float y, float z, double i)
    {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z),key)) return NULL;
      return addObservation(key,i);
    }

    // update inner nodes, sets color to average child color
    void updateInnerOccupancy();

    inline bool isNodeOccupied(const MeanOcTreeNode& node) const
    {
      return node.getMean() >= this->intensity_threshold_;
    }

    inline bool isNodeOccupied(const OcTreeNode* node) const
    {
      const MeanOcTreeNode* mean_node = reinterpret_cast<const MeanOcTreeNode*>(node);
      if (mean_node != NULL)
        return isNodeOccupied(*mean_node);
      else
        return false;
    }
    
  protected:

    double intensity_threshold_;

    void updateInnerOccupancyRecurs(MeanOcTreeNode* node, unsigned int depth);

    /**
     * Static member object which ensures that this OcTree's prototype
     * ends up in the classIDMapping only once. You need this as a 
     * static member in any derived octree class in order to read .ot
     * files through the AbstractOcTree factory. You should also call
     * ensureLinking() once from the constructor.
     */
    class StaticMemberInitializer{
       public:
         StaticMemberInitializer() {
           MeanOcTree* tree = new MeanOcTree(0.1);
           tree->clearKeyRays();
           AbstractOcTree::registerTreeType(tree);
         }

         /**
         * Dummy function to ensure that MSVC does not drop the
         * StaticMemberInitializer, causing this tree failing to register.
         * Needs to be called from the constructor of this octree.
         */
         void ensureLinking() {};
    };
    /// static member to ensure static initialization (only once)
    static StaticMemberInitializer meanOcTreeMemberInit;

  };

  //! user friendly output in format (r g b)
  std::ostream& operator<<(std::ostream& out, MeanOcTreeNode const* node);

} // end namespace

#endif
