import Vue from 'vue'

var version = "3";

var globalSearch = new Vue({
    el: '#global-search',
    data: {
        keyword: ''
    },
    methods: {
        onSearch: function () {
            if (globalSearch.keyword.length === 0) {
                location.href = '/';
                return;
            }
            location.href = '/?keyword=' + encodeURI(globalSearch.keyword);
        }
    }
});


const colorSeqView = new Vue({
    el: '#color-seq',
    data: {
        seq_id: '',
        sid: null,
        seqTexts: [],
        selectedIdToMotifs:{},
        motifTable: [],
        selectedSeq: {seq: ''},
        waitingBackend: true,
        picked: 'ss'
    },
    filters: {
        fixFloat: function (value) {
            value = Number(value);
            return value.toFixed(3);
        }
    },
    methods: {
        changeWaitingState: function(waitingState) {
            colorSeqView.waitingBackend = waitingState;
        },
        _updateMotifTable: function(row) {
            let motif = colorSeqView.motifTable[row];
            Vue.set(colorSeqView.motifTable, row, motif);
        },
        _fullUpdateMotifTable: function(seqInput) {
            colorSeqView.motifTable = [];

            let motifs = seqInput['motifs_16'];

            for (let i = 0; i < motifs.length; i++) {
                let offset = motifs[i]['offset'];

                let motifStr = '';
                if (i < motifs.length - 1) {
                    motifStr = seqInput['seq'].slice(offset, motifs[i+1]['offset']);
                }
                else {
                    motifStr = seqInput['seq'].slice(offset);
                }

                let motifRow = {
                    motifId: motifs[i]['id'],
                    row: i,
                    seq: motifStr,
                    offset: offset,
                    score: motifs[i]['score'],
                    tags: motifs[i]['tags'],
                    false_discovery: motifs[i]['false_discovery'],
                    manually_add: motifs[i]['manually_add']
                };
                Vue.set(colorSeqView.motifTable, i, motifRow);
            }
        },
        _updateSeqText: function(seqInput) {
            colorSeqView.seqTexts = [];

            let motifs = seqInput['motifs_16'];
            let colorPoses = new Set();
            for (let i = 0; i < motifs.length; ++i) {
                let motif = motifs[i];
                if (motif.false_discovery) {
                    continue;
                }
                for (let pos = motif.offset; pos < motif.offset + 16; ++pos) {
                    colorPoses.add(pos);
                }
            }

            let nsites = seqInput['nsites'];
            let nsitePoses = new Set();
            for (let i = 0; i < nsites.length; ++i) {
                let nsite = nsites[i];
                for (let pos = nsite.offset; pos < nsite.offset + 3; ++pos) {
                    nsitePoses.add(pos);
                }
            }

            let plainSeq = seqInput.seq;
            let ss = seqInput.ss.toLowerCase();
            for (let i = 0; i < plainSeq.length; ++i) {
                let seqText = {
                    offset: i,
                    seq: plainSeq[i],
                    ssStyle: 'ss' + ss[i],
                    innsite: nsitePoses.has(i),
                    inlrr: false
                };
                if (colorPoses.has(i)) {
                    seqText.inlrr = true;
                }
                Vue.set(colorSeqView.seqTexts, i, seqText);
            }
        },
        onSelectSeq: function (seqInput) {
            colorSeqView.selectedSeq = seqInput;
            colorSeqView.seq_id = seqInput['sequence_id'];
            colorSeqView.sid = seqInput['id'];
            colorSeqView._fullUpdateMotifTable(seqInput);
            colorSeqView._updateSeqText(seqInput);
        },
        markFalseDiscovery: function (row, falseDiscovery) {
            let motif = colorSeqView.motifTable[row];

            $.ajax({
                url: '/version/' + version + '/sequences/'+colorSeqView.sid+'/motifs/'+motif.motifId+'/false_discovery',
                method: 'PUT',
                contentType: "application/json",
                data: JSON.stringify({false_discovery: falseDiscovery}),
                success: function (response) {
                    motif.false_discovery = falseDiscovery;
                    colorSeqView._updateMotifTable(row);
                    colorSeqView._updateSeqText(colorSeqView.selectedSeq);
                },
                error: function (xhr, textStatus, errorThrown) {
                    if (xhr.status === 429) {
                        $('#tooManyRequestAlert').modal('show');
                        return;
                    }
                    let errorBody = JSON.parse(xhr.responseText);
                    alert(errorThrown + ': ' + errorBody['message']);
                }
            });
        },
        removeMotif: function(motifId) {
            $.ajax({
                url: '/version/' + version + '/sequences/'+colorSeqView.sid+'/motifs/'+motifId,
                type: 'DELETE',
                contentType: "application/json",
                dataType: "json",
                success: function (response) {
                    for (let i = 0; i < seqsRow.activeSeq.motifs_16.length; ++i) {
                        if (seqsRow.activeSeq.motifs_16[i].id === motifId) {
                            seqsRow.activeSeq.motifs_16.splice(i, 1);
                            break;
                        }
                    }
                    colorSeqView.onSelectSeq(seqsRow.activeSeq);
                },
                error: function (xhr, textStatus, errorThrown) {
                    let errorBody = JSON.parse(xhr.responseText);
                    alert(errorThrown + ': ' + errorBody['message']);
                }
            });

        },
        markLrr: function () {
            if (window.getSelection().anchorNode === null) {
                alert("Please select the sequence first, and the offset of the sequence will be marked as the beginning of the LRR. The length of the LRR is 16, regardless of the selected length.");
                return;
            }
            let offset = window.getSelection().anchorNode.parentElement.getAttribute('offset');
            if (offset === null) {
                alert("It is currently not supported to select a sequence in the table. Please select it in the sequence below.");
                return;
            }

            $.ajax({
                url: '/version/' + version + '/sequences/'+colorSeqView.sid+'/motifs',
                method: 'POST',
                contentType: "application/json",
                data: JSON.stringify({offset: parseInt(offset)}),
                success: function (response) {
                    let motif = response.motifs_16[0];
                    console.log(motif);
                    seqsRow.activeSeq.motifs_16.push({
                        id: motif.id,
                        offset: motif.offset,
                        score: motif.score,
                        tags: new Set(motif.tags),
                        manually_add: motif.manually_add
                    });
                    seqsRow.activeSeq.motifs_16.sort(function(a,b){return a.offset - b.offset});
                    colorSeqView.onSelectSeq(seqsRow.activeSeq);
                },
                error: function (xhr, textStatus, errorThrown) {
                    let errorBody = JSON.parse(xhr.responseText);
                    alert(errorThrown + ': ' + errorBody['message']);
                }
            });
        }
    }
});

const seqsRow = new Vue({
    el: '#navigator-view',
    components: {
        'multiselect': window.VueMultiselect.default
    },
    data: {
        total: 0,
        activeSeq: null,
        seqs: [],
        keyword: "",
        pageIndex: -1,
        pageSize: -1,
        goToIndex: 0,
        filterNamesToValue: {},
        waitingBackend: false,
        value: null,
        options: [
            {code: 'AMBTR', name: 'Amborella trichopoda'},
            {code: 'ARALY', name: 'Arabidopsis lyrata'},
            {code: 'ARATH', name: 'Arabidopsis thaliana'},
            {code: 'BRADI', name: 'Brachypodium distachyon'},
            {code: 'BRARA', name: 'Brassica rapa'},
            {code: 'GLYMA', name: 'Glycine max '},
            {code: 'MARPL', name: 'Marchantia polymorpha'},
            {code: 'MEDTR', name: 'Medicago truncatula'},
            {code: 'ORYSI', name: 'Oryza sativa ssp. Indica'},
            {code: 'ORYSJ', name: 'Oryza sativa ssp. Japonica'},
            {code: 'PHODC', name: 'Phoenix dactylifera'},
            {code: 'PHYPA', name: 'Physcomitrella patens'},
            {code: 'POPTR', name: 'Populus trichocarpa'},
            {code: 'SELML', name: 'Selaginella moellendorffii'},
            {code: 'SOLLC', name: 'Solanum lycopersicum'},
            {code: 'SOLTU', name: 'Solanum tuberosum'},
            {code: 'MAIZE', name: 'Zea mays'},
        ]
    },
    methods: {
        _selectSeq: function(seq) {
            seq['isActive'] = true;
            colorSeqView.onSelectSeq(seq)
        },
        onSelect: function (seq) {
            if (seqsRow.activeSeq === seq) {
                return;
            }
            seqsRow.activeSeq['isActive'] = false;
            seqsRow.activeSeq = seq;
            seqsRow._selectSeq(seq);
        },
        _changeWaitingState: function(isWating) {
            seqsRow.waitingBackend = isWating;
            colorSeqView.changeWaitingState(isWating);
        },
        changePage: function (page, size, filterNamesToValue) {
            seqsRow.pageIndex = page;
            seqsRow.goToIndex = seqsRow.pageIndex + 1;
            seqsRow.pageSize = size;
            seqsRow.filterNamesToValue = filterNamesToValue;
            seqsRow._changeWaitingState(true);

            let searchParas = $.extend({}, filterNamesToValue);
            searchParas.page = page;
            searchParas.size = size;
            seqsRow.keyword = seqsRow.keyword.split(' ').join('');
            if (seqsRow.keyword.length > 0) {
                searchParas['keyword'] = seqsRow.keyword;
            }

            $.ajax({
                url: '/version/' + version + '/sequences',
                method: 'GET',
                contentType: "application/json",
                data: searchParas,
                success: function (data) {
                    seqsRow._changeWaitingState(false);
                    seqsRow.seqs = [];
                    let first = true;
                    for (let i = 0; i < data.sequences.length; i++) {
                        let seq = data.sequences[i];
                        let motifs = seq['motifs_16'];
                        for (let i = 0; i < motifs.length; i++) {
                            motifs[i]['tags'] = new Set(motifs[i]['tags']);
                        }
                        if (i === 0) {
                            seqsRow.activeSeq = seq;
                            seqsRow._selectSeq(seq);
                            first = false;
                        } else {
                            seq['isActive'] = false;
                        }
                        Vue.set(seqsRow.seqs, i, seq);
                    }
                    seqsRow.total = data.total;
                },
                error: function (xhr, textStatus, errorThrown) {
                    seqsRow._changeWaitingState(false);
                    if (xhr.status === 429) {
                        $('#tooManyRequestAlert').modal('show');
                    }
                }
            });
        },
        onSearch: function() {
            if (seqsRow.keyword.length === 0) {
                return;
            }
            seqsRow.changePage(0, seqsRow.pageSize, seqsRow.filterNamesToValue);
        },
        onSpeciesSelectEnd(value, id) {
            if (seqsRow.value === null) {
                delete seqsRow.filterNamesToValue.species;
            } else {
                seqsRow.filterNamesToValue['species'] = [seqsRow.value.code];
            }
            seqsRow.changePage(0, seqsRow.pageSize, seqsRow.filterNamesToValue);

        },
        prevPage: function () {
            if (seqsRow.pageIndex > 0) {
                seqsRow.changePage(seqsRow.pageIndex - 1, seqsRow.pageSize, seqsRow.filterNamesToValue);
            }
        },
        nextPage: function () {
            if (seqsRow.pageIndex < Math.ceil(seqsRow.total/seqsRow.pageSize) - 1) {
                seqsRow.changePage(seqsRow.pageIndex + 1, seqsRow.pageSize, seqsRow.filterNamesToValue);
            }
        },
        goTo: function() {
            let index = seqsRow.goToIndex - 1;
            if (index >= 0 && index <= Math.ceil(seqsRow.total/seqsRow.pageSize) - 1) {
                seqsRow.changePage(index, seqsRow.pageSize, seqsRow.filterNamesToValue);
            }
        },
        onChangePageSize: function () {
            seqsRow.changePage(seqsRow.pageIndex, seqsRow.pageSize, seqsRow.filterNamesToValue);
        }
    }
});


const filtersArea = new Vue({
    el: '#filters-area',
    data: {
        offsetlt: null,
        offsetgt: null,
        lrr_count_lt: null,
        lrr_count_gt: null,
        keyword: null
    },
    methods: {
        _get_offset: function(filters) {
            if (filtersArea.offsetlt != null
                && filtersArea.offsetgt
                && filtersArea.offsetlt === filtersArea.offsetgt) {
                filters['offseteq'] = filtersArea.offsetlt;
                return;
            }
            if (filtersArea.offsetlt != null){
                filters.offsetlt = filtersArea.offsetlt;
            }
            if (filtersArea.offsetgt != null) {
                filters['offsetgt'] = filtersArea.offsetgt;
            }
        },
        _get_lrr_count: function(filters) {
            if (filtersArea.lrr_count_lt != null
                && filtersArea.lrr_count_gt
                && filtersArea.lrr_count_lt === filtersArea.lrr_count_gt) {
                filters['lrr_count_eq'] = filtersArea.lrr_count_lt;
                return;
            }
            if (filtersArea.lrr_count_lt != null){
                filters['lrr_count_lt'] = filtersArea.lrr_count_lt;
            }
            if (filtersArea.lrr_count_gt != null) {
                filters['lrr_count_gt'] = filtersArea.lrr_count_gt;
            }
        },
        onSearch: function() {
            let filters = {};
            filtersArea._get_offset(filters);
            filtersArea._get_lrr_count(filters);
            if (filtersArea.keyword != null){
                filters['keyword'] = filtersArea.keyword;
            }
            seqsRow.changePage(0, seqsRow.pageSize, filters);
        }
    }
});


const getAllQueries = function() {
    let searchString = window.location.search.slice(1);
    let searchPairs = searchString.split('&');
    let filterNamesToValue = {};
    for (let i = 0; i < searchPairs.length; ++i) {
        let kv = searchPairs[i].split('=');
        if (kv[0].length === 0) {
            continue;
        }
        if (kv[1] === undefined || kv[1].length === 0) {
            filterNamesToValue[kv[0]] = true;
        }
        else {
            filterNamesToValue[kv[0]] = decodeURIComponent(kv[1]);
        }
    }
    return filterNamesToValue;
};

let path = window.location.pathname.toLowerCase();
if (path.substring(0, 3) === "/v1") {
    version = "1";
} else if (path.substring(0, 3) === "/v2") {
    version = "2";
}

seqsRow.changePage(0, 20, getAllQueries());