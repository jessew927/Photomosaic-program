/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
  if(first[curDim] < second[curDim])
    return true;
  else if(first[curDim] > second[curDim])
    return false;
  else
    return first < second;

}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
  double curBestDistance = 0;
  double potentialDistance = 0;

  for(int i = 0; i < Dim; i++){
    curBestDistance = curBestDistance + pow((target[i] - currentBest[i]), 2);
    potentialDistance = potentialDistance + pow((target[i] - potential[i]), 2);
  }

  if(potentialDistance < curBestDistance)
    return true;
  else if(potentialDistance > curBestDistance)
    return false;
  else
    return potential < currentBest;
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
  if(newPoints.empty()){
    size = 0;
    root = nullptr;
  }
  else{
    size = newPoints.size();
    vector<Point<Dim>> copyNewPoints;
    copyNewPoints.assign(newPoints.begin(), newPoints.end());
    root = constructorHelper(copyNewPoints, 0, 0, newPoints.size()-1);
  }
}

template <int Dim>
typename KDTree<Dim>::KDTreeNode* KDTree<Dim>::constructorHelper(vector<Point<Dim>>& v, int curDim, int left, int right){
  if(left == right){
    KDTreeNode * node = new KDTreeNode(v[left]);
    return node;
  }
  else{
    int medianIndex = ceil((left+right)/2);
    Point<Dim> medianPoint = quickSelect(v, left, right, medianIndex, curDim);
    KDTreeNode * subRoot = new KDTreeNode(medianPoint);
    curDim = (curDim + 1) % Dim;
    if(left < medianIndex)
      subRoot->left = constructorHelper(v, curDim, left, medianIndex - 1);
    if(right > medianIndex)
      subRoot->right = constructorHelper(v, curDim, medianIndex + 1, right);
    return subRoot;
  }
}

template<int Dim>
Point<Dim> KDTree<Dim>::quickSelect(vector<Point<Dim>>& v, int left, int right, int k, int curDim){
  if(left == right)
    return v[left];

  int pivot = partitionHelper(v, left, right, curDim);

  if(k == pivot)
    return v[k];
  else if(k < pivot)
    return quickSelect(v, left, pivot - 1, k, curDim);
  else
    return quickSelect(v, pivot + 1, right, k , curDim);
}

template<int Dim>
int KDTree<Dim>::partitionHelper(vector<Point<Dim>> & v, int left, int right, int curDim){
  Point<Dim> rightPoint = v[right];
  int leftIndex = left;
  for(int i = left; i <= right; i++){
    if(smallerDimVal(v[i], rightPoint, curDim)){
      swap(v[i], v[leftIndex]);
      leftIndex++;
    }
  }

  swap(v[leftIndex], v[right]);
  return leftIndex;
}

template<int Dim>
void KDTree<Dim>::swap(Point<Dim>& x, Point<Dim>& y){
  Point<Dim> temp = x;
  x = y;
  y = temp;
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  copy(other);
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  if(this != &rhs){
    destroy(this->root);
    copy(rhs);
  }

  return *this;
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  destroy(this->root);
}

template<int Dim>
void KDTree<Dim>::copy(const KDTree<Dim>& copyTree){
  size = copyTree.size;
  root = new KDTreeNode(copyTree.root->point);
  copyTreeNodes(root, copyTree.root);
}

template<int Dim>
void KDTree<Dim>::copyTreeNodes(KDTreeNode* curNode, KDTreeNode* other){
  if(other == nullptr)
    return;
  else{
    if(other->left != nullptr)
      curNode->left = new KDTreeNode(other->left->point);
    if(other->right != nullptr)
      curNode->right = new KDTreeNode(other->right->point);

    copyTreeNodes(curNode->left, other->left);
    copyTreeNodes(curNode->right, other->right);
  }
}

template<int Dim>
void KDTree<Dim>::destroy(KDTreeNode * curNode){
  if(curNode == nullptr)
    return;
  else{
    destroy(curNode->left);
    destroy(curNode->right);
    delete curNode;
  }
}


template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
  KDTreeNode * curNode = this->root;
  Point<Dim> nearestPoint;
  double bestRadius = 0;
  bool firstRun = true;
  FNNHelper(curNode, query, nearestPoint, 0, bestRadius, firstRun);
  return nearestPoint;
}

template <int Dim>
void KDTree<Dim>::FNNHelper(KDTreeNode* curNode, const Point<Dim>& query ,Point<Dim>& currentBest, int curDim, double& bestRadius, bool& firstRun) const{
  bool visitedLeft = false, vistedRight = false;
  int nextDim = (curDim + 1) % Dim;
  if(curNode == nullptr)
    return;
  if(curNode->left == nullptr && curNode->right == nullptr){
    if(firstRun){
      currentBest = curNode->point;
      firstRun = false;
    }
    else{
      if(shouldReplace(query, currentBest, curNode->point))
        currentBest = curNode->point;
    }
    bestRadius = calculateDistance(query, currentBest);
    return;
  }
  else if(curNode->left == nullptr && curNode->right != nullptr){
    FNNHelper(curNode->right, query, currentBest, nextDim, bestRadius, firstRun);
    vistedRight = true;
  }
  else if(curNode->left != nullptr && curNode->right == nullptr){
    FNNHelper(curNode->left, query, currentBest, nextDim, bestRadius, firstRun);
    visitedLeft = true;
  }
  else{
    if(curNode->point == query){
       currentBest = curNode->point;
     }
    else if(smallerDimVal(query,curNode->point, curDim)){
      FNNHelper(curNode->left, query, currentBest, nextDim, bestRadius, firstRun);
      visitedLeft = true;
    }
    else if(!smallerDimVal(query,curNode->point, curDim)){
      FNNHelper(curNode->right, query, currentBest, nextDim, bestRadius, firstRun);
      vistedRight = true;
    }
  }

  int potentialDistance = calculateDistance(query, curNode->point);
  if(shouldReplace(query, currentBest, curNode->point)){
    currentBest = curNode->point;
    bestRadius = potentialDistance;
  }

  int splitDistance = pow(curNode->point[curDim] - query[curDim], 2);

  if(splitDistance <= bestRadius){
    if(visitedLeft)
      FNNHelper(curNode->right, query, currentBest, nextDim, bestRadius, firstRun);
    else if(vistedRight)
      FNNHelper(curNode->left, query, currentBest, nextDim, bestRadius, firstRun);
  }

}

template <int Dim>
double KDTree<Dim>::calculateDistance(const Point<Dim>& query ,Point<Dim>& curPoint) const{
  double curDistance = 0;
  for(int i = 0; i < Dim; i++){
    curDistance = curDistance + pow((query[i] - curPoint[i]), 2);
  }
  return curDistance;
}
