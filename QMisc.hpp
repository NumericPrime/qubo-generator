#include<unordered_map>
#include<utility>
size_t QMiscCantorMapping(size_t x,size_t y){
    return (((x+y+1)*(x+y))>>1)+y;
}
template<> struct std::hash<std::pair<unsigned,unsigned>>{
    std::size_t operator()(const std::pair<unsigned,unsigned> &f) const{
        return QMiscCantorMapping(f.first,f.second);
    }
};