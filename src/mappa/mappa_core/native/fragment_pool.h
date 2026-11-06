// A reaction's unique completed fragments certify smaller existence queries.
// No growth trajectory, deferred boundary or branch outcome is cached here.
#include "fragment_cache_io.h"

class FragmentQuery;
class FragmentPool {
    struct Entry;
    using Posting = std::list<Entry*>;
    struct Entry {
        std::string key;
        std::shared_ptr<const Fragment> fragment;
        std::vector<int> image;
        std::vector<Posting::iterator> positions;
        size_t bytes=0;
        bool imported=false;
    };
    using Entries=std::list<Entry>;
    Entries entries;
    std::unordered_map<std::string,Entries::iterator> index;
    std::vector<Posting> by_atom;
    size_t budget, resident=0, peak=0;
    long lookups=0,hits=0,checks=0,evictions=0,oversized=0;
    uint64_t revision=0;

    template<class T> static void append(std::string& key,const T& value) {
        key.append(reinterpret_cast<const char*>(&value),sizeof(value));
    }
    static std::string key(const Fragment& fragment) {
        std::string result;
        for (const auto* atoms:{&fragment.source->atoms,&fragment.target->atoms}) {
            append(result,atoms->size());
            for (int atom:*atoms) append(result,atom);
        }
        return result;
    }
    static size_t bytes(const Fragment& fragment) {
        // Conservatively charge shared self records to each owning fragment.
        size_t result=sizeof(Fragment)+sizeof(Entry)+512;
        for (const auto* self:{fragment.source.get(),fragment.target.get()}) {
            result+=sizeof(PreparedSymmetry)+self->weights.capacity()*sizeof(double)
                +(self->atoms.capacity()+self->anchors.capacity())*sizeof(int)
                +self->elements.capacity()*sizeof(std::string)
                +self->generators.capacity()*sizeof(std::vector<int>);
            for (const auto& element:self->elements) result+=element.capacity()+1;
            for (const auto& generator:self->generators) result+=generator.capacity()*sizeof(int);
        }
        result+=(fragment.anchors.capacity()+fragment.accepted_edges.capacity()+fragment.target_edges.capacity())*sizeof(Pair)
            +fragment.witness.capacity()*sizeof(int)+fragment.edge_tolerances.capacity()*sizeof(double)
            +fragment.domains.capacity()*sizeof(std::vector<int>);
        for (const auto& domain:fragment.domains) result+=domain.capacity()*sizeof(int);
        return result;
    }
    static bool same_record(const Fragment& a,const Fragment& b) {
        return a.witness==b.witness && a.domains==b.domains && a.anchors==b.anchors
            && a.accepted_edges==b.accepted_edges && a.edge_tolerances==b.edge_tolerances
            && a.target_edges==b.target_edges && a.inclusive==b.inclusive
            && a.target_bond_floor==b.target_bond_floor
            && a.source->bond_threshold==b.source->bond_threshold
            && a.target->bond_threshold==b.target->bond_threshold
            && a.source->anchors==b.source->anchors && a.target->anchors==b.target->anchors;
    }
    void remove_oldest() {
        auto& entry=entries.front();
        for (size_t i=0;i<entry.fragment->source->atoms.size();++i)
            by_atom[entry.fragment->source->atoms[i]].erase(entry.positions[i]);
        resident-=entry.bytes;index.erase(entry.key);entries.pop_front();++evictions;
    }
public:
    Source base,active_source;
    const Target& target;
    std::vector<Pair> active_cuts;
    FragmentPool(const Source& source,const Target& target,size_t budget)
        :by_atom(source.n),budget(budget),base(source),active_source(source),target(target) {}

    std::shared_ptr<const Fragment> remember(std::shared_ptr<const Fragment> fragment) {
        auto identity=key(*fragment);
        auto found=index.find(identity);
        if (found!=index.end()) {
            entries.splice(entries.end(),entries,found->second);
            // Share equal scientific records. A new decision keeps its own
            // locks and accepted edges even when its atom-set identity repeats.
            return same_record(*found->second->fragment,*fragment) ? found->second->fragment : fragment;
        }
        size_t charge=bytes(*fragment)+2*(identity.capacity()+1)+base.n*sizeof(int)
            +fragment->source->atoms.size()*(sizeof(Posting::iterator)+64);
        if (charge>budget) {++oversized;return fragment;}
        while (!entries.empty() && resident+charge>budget) remove_oldest();
        entries.emplace_back();auto entry=std::prev(entries.end());
        entry->key=std::move(identity);entry->fragment=fragment;entry->bytes=charge;
        entry->image.assign(base.n,-1);
        for (size_t i=0;i<fragment->source->atoms.size();++i) {
            int atom=fragment->source->atoms[i];
            entry->image[atom]=fragment->target->atoms[fragment->witness[i]];
            by_atom[atom].push_back(&*entry);
            entry->positions.push_back(std::prev(by_atom[atom].end()));
        }
        index.emplace(entry->key,entry);resident+=charge;peak=std::max(peak,resident);++revision;
        return fragment;
    }

    std::vector<int> find(const FragmentQuery& query);

    void cut(std::vector<Pair> cuts) {
        for (auto& [a,b]:cuts) if (a>b) std::swap(a,b);
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        if (cuts==active_cuts) return;
        for (auto [a,b]:cuts)
            if (a<0 || b<0 || a>=base.n || b>=base.n || !base.has_edge(a,b))
                throw std::invalid_argument("cut is not an active edge of the source");
        active_source=base;active_cuts=std::move(cuts);
        for (auto [a,b]:active_cuts) {
            active_source.edge[(size_t)a*base.n+b]=active_source.edge[(size_t)b*base.n+a]=0;
            for (auto [x,y]:{Pair{a,b},Pair{b,a}}) {
                auto& neighbors=active_source.adj[x];
                neighbors.erase(std::remove(neighbors.begin(),neighbors.end(),y),neighbors.end());
            }
        }
    }

    py::bytes snapshot_context() const {
        fragment_cache_io::Writer out;
        out(fragment_cache_io::schema,uint64_t(sizeof(long)));
        for (const Graph* graph:{static_cast<const Graph*>(&base),static_cast<const Graph*>(&target)})
            out(graph->n,graph->elem,graph->adj,graph->edge,graph->w,graph->bond_cut,
                graph->bond_tolerances,graph->tolerance_inclusive);
        return py::bytes(out.data);
    }
    py::bytes dump(const std::string& identity) const {
        fragment_cache_io::Writer out;
        out(fragment_cache_io::schema,identity,uint64_t(entries.size()));
        for (const auto& entry:entries) out(entry.fragment);
        return py::bytes(out.data);
    }
    bool load(const std::string& data,const std::string& identity) {
        fragment_cache_io::Reader in(data);
        std::string schema,observed;in(schema,observed);
        if (schema!=fragment_cache_io::schema) throw std::invalid_argument("Unsupported fragment pool schema");
        if (observed!=identity) throw std::invalid_argument("Fragment pool context differs");
        std::vector<std::shared_ptr<const Fragment>> values(in.size());
        for (auto& value:values) in(value);
        in.finish();
        for (const auto& value:values) {
            if (value->source->atoms.size()!=value->witness.size())
                throw std::invalid_argument("Fragment witness size differs");
            std::vector<uint8_t> used(target.n,0);
            for (size_t i=0;i<value->witness.size();++i) {
                int r=value->source->atoms[i],local=value->witness[i];
                if (r<0 || r>=base.n || local<0 || local>=int(value->target->atoms.size()))
                    throw std::invalid_argument("Fragment witness index outside endpoints");
                int p=value->target->atoms[local];
                if (p<0 || p>=target.n || used[p] || base.ecode[r]!=target.ecode[p])
                    throw std::invalid_argument("Fragment witness is not an injective element match");
                used[p]=1;
            }
        }
        for (auto& entry:entries) entry.imported=false;
        for (const auto& value:values) {
            remember(value);
            auto found=index.find(key(*value));
            if (found!=index.end()) found->second->imported=true;
        }
        return std::any_of(entries.begin(),entries.end(),[](const Entry& entry){return !entry.imported;});
    }
    py::dict stats() const {
        py::dict out;
        out["existence_queries"]=lookups;out["fragment_hits"]=hits;out["witness_checks"]=checks;
        out["entries"]=entries.size();out["resident_bytes"]=resident;out["peak_resident_bytes"]=peak;
        out["cache_budget_bytes"]=budget;out["evictions"]=evictions;out["oversized_results"]=oversized;
        out["revision"]=revision;
        return out;
    }
};
