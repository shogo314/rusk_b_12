template<class T>
struct RuskPersistentSegTreeNode {
    Obj<T> value;
    Obj<RuskPersistentSegTreeNode> left;
    Obj<RuskPersistentSegTreeNode> right;
};
