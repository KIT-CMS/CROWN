namespace metfilter{

    auto ApplyMetFilter(auto df, std::string filtername){
        return df.Filter([](const bool flag){return flag;},{filtername});
    }

} // end namespace metfilters