// A reaction's unique completed fragments certify smaller existence queries.
// No growth trajectory, deferred boundary or branch outcome is cached here.
#include "fragment_cache_io.h"

class FragmentPool {
    struct Entry {
        std::string key;
        std::shared_ptr<const Fragment> fragment;
        std::vector<int> image;
        size_t block=0;
        int slot=0;
        size_t bytes=0;
        bool imported=false;
    };
    // Each bit names one stored witness. Intersections test 64 fragments at once.
    struct Block {
        uint64_t live=0;
        Entry* entries[64]{};
        std::vector<uint64_t> atoms, edges, assignments;
        Block(size_t source, size_t target, size_t bonds)
            :atoms(source),edges(bonds),assignments(source*target) {}
    };
    using Entries=std::list<Entry>;
    Entries entries;
    std::unordered_map<std::string,Entries::iterator> index;
    std::vector<std::unique_ptr<Block>> blocks;
    std::vector<Pair> bonds;
    std::vector<int> bond_index;
    std::vector<size_t> atom_counts,edge_counts;
    std::vector<uint64_t> atom_blocks,edge_blocks;
    size_t block_words=0,index_bytes=0;
    std::vector<int> query_atoms,query_edges,query_outside;
    std::optional<double> indexed_tolerance;
    size_t budget, resident=0, peak=0;
    long lookups=0,hits=0,block_checks=0,evictions=0,oversized=0;
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
    size_t block_bytes() const {
        return sizeof(Block)+128+sizeof(uint64_t)*(base.n+bonds.size()+size_t(base.n)*target.n);
    }
    size_t available_block() const {
        for (size_t i=0;i<blocks.size();++i)
            if (!blocks[i] || blocks[i]->live!=~uint64_t(0)) return i;
        return blocks.size();
    }
    // Charge vector growth capacity as well as occupied index words.
    size_t index_word_bytes() const {return 2*(size_t(base.n)+bonds.size())*sizeof(uint64_t);}
    void grow_index() {
        ++block_words;
        atom_blocks.resize(size_t(base.n)*block_words);
        edge_blocks.resize(bonds.size()*block_words);
        resident+=index_word_bytes();index_bytes+=index_word_bytes();
    }
    void index_edges(const Entry& entry) {
        if (!indexed_tolerance) return;
        auto& block=*blocks[entry.block];
        uint64_t bit=uint64_t(1)<<entry.slot;
        for (size_t i=0;i<bonds.size();++i) {
            auto [a,b]=bonds[i];
            int p=entry.image[a],q=entry.image[b];
            if (p>=0 && q>=0 && target.has_edge(p,q) && growth_edge_supported(
                    base.wbo(a,b),target.wbo(p,q),base.tolerance(a,b,*indexed_tolerance),
                    target.bond_cut,base.tolerance_inclusive)) {
                block.edges[i]|=bit;++edge_counts[i];
                edge_blocks[(entry.block/64)*bonds.size()+i]|=uint64_t(1)<<(entry.block%64);
            }
        }
    }
    void prepare_tolerance(double tolerance) {
        if (indexed_tolerance && *indexed_tolerance==tolerance) return;
        indexed_tolerance=tolerance;
        std::fill(edge_counts.begin(),edge_counts.end(),0);
        std::fill(edge_blocks.begin(),edge_blocks.end(),0);
        for (auto& block:blocks) if (block) std::fill(block->edges.begin(),block->edges.end(),0);
        for (const auto& entry:entries) index_edges(entry);
    }
    void remove_oldest() {
        auto& entry=entries.front();
        auto& block=*blocks[entry.block];
        uint64_t keep=~(uint64_t(1)<<entry.slot);
        for (int r:entry.fragment->source->atoms) {
            block.atoms[r]&=keep;
            --atom_counts[r];
            if (!block.atoms[r]) atom_blocks[(entry.block/64)*base.n+r]&=~(uint64_t(1)<<(entry.block%64));
            block.assignments[size_t(r)*target.n+entry.image[r]]&=keep;
        }
        for (size_t i=0;i<bonds.size();++i) {
            if (block.edges[i]&~keep) --edge_counts[i];
            block.edges[i]&=keep;
            if (!block.edges[i]) edge_blocks[(entry.block/64)*bonds.size()+i]&=~(uint64_t(1)<<(entry.block%64));
        }
        block.entries[entry.slot]=nullptr;
        block.live&=keep;
        if (!block.live) {blocks[entry.block].reset();resident-=block_bytes();}
        resident-=entry.bytes;index.erase(entry.key);entries.pop_front();++evictions;
    }
public:
    Source base,active_source;
    const Target& target;
    std::vector<Pair> active_cuts;
    FragmentPool(const Source& source,const Target& target,size_t budget)
        :bond_index(size_t(source.n)*source.n,-1),atom_counts(source.n),budget(budget),
         base(source),active_source(source),target(target) {
        for (int a=0;a<base.n;++a) for (int b:base.adj[a]) if (a<b) {
            bond_index[size_t(a)*base.n+b]=bond_index[size_t(b)*base.n+a]=int(bonds.size());
            bonds.push_back({a,b});
        }
        edge_counts.resize(bonds.size());
    }

    std::shared_ptr<const Fragment> remember(std::shared_ptr<const Fragment> fragment) {
        auto identity=key(*fragment);
        auto found=index.find(identity);
        if (found!=index.end()) {
            entries.splice(entries.end(),entries,found->second);
            // Share equal scientific records. A new decision keeps its own
            // locks and accepted edges even when its atom-set identity repeats.
            return same_record(*found->second->fragment,*fragment) ? found->second->fragment : fragment;
        }
        size_t charge=bytes(*fragment)+2*(identity.capacity()+1)+base.n*sizeof(int);
        if (charge+block_bytes()+std::max(index_bytes,index_word_bytes())>budget) {++oversized;return fragment;}
        size_t number;
        while (true) {
            number=available_block();
            size_t extra=number==blocks.size() || !blocks[number] ? block_bytes() : 0;
            if (number/64>=block_words) extra+=index_word_bytes();
            if (resident+charge+extra<=budget) break;
            remove_oldest();
        }
        if (number/64>=block_words) grow_index();
        if (number==blocks.size()) blocks.emplace_back();
        if (!blocks[number]) {
            blocks[number]=std::make_unique<Block>(base.n,target.n,bonds.size());
            resident+=block_bytes();
        }
        auto& block=*blocks[number];
        entries.emplace_back();auto entry=std::prev(entries.end());
        entry->key=std::move(identity);entry->fragment=fragment;entry->bytes=charge;
        entry->block=number;entry->slot=__builtin_ctzll(~block.live);
        uint64_t bit=uint64_t(1)<<entry->slot;
        block.live|=bit;block.entries[entry->slot]=&*entry;
        entry->image.assign(base.n,-1);
        for (size_t i=0;i<fragment->source->atoms.size();++i) {
            int atom=fragment->source->atoms[i];
            entry->image[atom]=fragment->target->atoms[fragment->witness[i]];
            block.atoms[atom]|=bit;
            ++atom_counts[atom];
            atom_blocks[(number/64)*base.n+atom]|=uint64_t(1)<<(number%64);
            block.assignments[size_t(atom)*target.n+entry->image[atom]]|=bit;
        }
        index_edges(*entry);
        index.emplace(entry->key,entry);resident+=charge;peak=std::max(peak,resident);++revision;
        return fragment;
    }

    std::vector<int> find(const std::vector<uint8_t>& selected,const std::vector<Pair>& edges,
                          const std::vector<int>& locked,double tolerance) {
        ++lookups;
        if (entries.empty()) return {};
        prepare_tolerance(tolerance);
        query_atoms.clear();query_edges.clear();query_outside.clear();
        size_t least=entries.size()+1;
        int rare_atom=-1,rare_edge=-1;
        for (auto [a,b]:edges) {
            int edge=bond_index[size_t(a)*base.n+b];
            if (edge<0 || !edge_counts[edge]) return {};
            query_edges.push_back(edge);
            if (edge_counts[edge]<least) {least=edge_counts[edge];rare_edge=edge;}
        }
        for (int r=0;r<base.n;++r) {
            if (selected[r]) {
                if (!atom_counts[r]) return {};
                query_atoms.push_back(r);
                if (atom_counts[r]<least) {least=atom_counts[r];rare_atom=r;rare_edge=-1;}
            } else if (locked[r]>=0) query_outside.push_back(locked[r]);
        }
        if (query_atoms.empty()) return {};
        const auto& active=rare_edge>=0 ? edge_blocks : atom_blocks;
        size_t feature=rare_edge>=0 ? rare_edge : rare_atom;
        size_t width=rare_edge>=0 ? bonds.size() : size_t(base.n);
        for (size_t word=0;word<(blocks.size()+63)/64;++word) {
            uint64_t pending=active[word*width+feature];
            while (pending) {
                size_t number=64*word+__builtin_ctzll(pending);pending&=pending-1;
                const auto& block=blocks[number];
                ++block_checks;
                uint64_t candidates=rare_edge>=0 ? block->edges[rare_edge] : block->atoms[rare_atom];
                if (!candidates) continue;
                for (int edge:query_edges) {candidates&=block->edges[edge];if (!candidates) break;}
                if (!candidates) continue;
                for (int r:query_atoms) {candidates&=block->atoms[r];if (!candidates) break;}
                if (!candidates) continue;
                // Anchors intersect assignment sets; occupied targets exclude them.
                // Injectivity of each stored witness handles locks inside the query.
                for (int r:query_atoms) if (locked[r]>=0) {
                    candidates&=block->assignments[size_t(r)*target.n+locked[r]];
                    if (!candidates) break;
                }
                if (!candidates) continue;
                for (int r:query_atoms) {
                    for (int p:query_outside) candidates&=~block->assignments[size_t(r)*target.n+p];
                    if (!candidates) break;
                }
                if (!candidates) continue;
                auto* entry=block->entries[__builtin_ctzll(candidates)];
                ++hits;
                entries.splice(entries.end(),entries,index.at(entry->key));
                std::vector<int> witness(base.n,-1);
                for (int r:query_atoms) witness[r]=entry->image[r];
                return witness;
            }
        }
        return {};
    }

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
        out["existence_queries"]=lookups;out["fragment_hits"]=hits;out["set_blocks_checked"]=block_checks;
        out["entries"]=entries.size();out["resident_bytes"]=resident;out["peak_resident_bytes"]=peak;
        out["cache_budget_bytes"]=budget;out["evictions"]=evictions;out["oversized_results"]=oversized;
        out["revision"]=revision;
        return out;
    }
};
