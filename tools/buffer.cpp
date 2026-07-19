void ofBuffer::append(const string &buffer);
void ofBuffer::RLines::getReverseLines();
class Init_PWM: public Motif{
    private:
        score_mat mx;
        std::vector<unsigned int> idx;
        std::vector<double> score;
        unsigned int len;
        unsigned int win;
        unsigned int a_size;
        unsigned int win_idx;
        double thres;
    public:
        Init_PWM(const score_mat& mx, const vector<double>& vec, unsigned int win, double thres);
        std::pair<bool, double> win_match()

};
class Aggregate_PWM: public Motif {
    private:

}
void score_mat reverse_complement(const score_mat &mx){
    size_t a=mx.size();
    size_t n=mx[0].size();

}
