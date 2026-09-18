
const TogoCoord = () => {
  let config = false;
  let id2label = false;
  const api = "https://sparql-support.dbcls.jp/sparqlist/api/";
  let options = {
    method: 'POST',
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/x-www-form-urlencoded'
    }
  }
  const checkMode = () => {
    let id = false;
    document.querySelectorAll(".tab").forEach(element => {
      if (element.classList.contains("tab_selected")) {
	id = element.id;
      }
    });
    return id;
  }

  const checkQueryCategory = (query, config) => {
    let new_query = query;
    let category = false;
    if (query.match(/^[a-z_]+#/)) {
      const string = query.match(/^([a-z_]+)#/)[1];
      for (let d of config.categories) {
	if (d.id == string) {
	  category = d.id;
	  new_query = query.match(/^[a-z_]+#(.+)$/)[1];
	  break;
	}
      }
    }
    return [new_query, category];
  }

  const getQueryCategories = (query) => {
    return ["transcript", "mrna", "protein"];
  }
  
  const showCategoryList = (depth, category, query, categories) => {
    //console.log(query);
    const output = document.querySelector("#nested_div_" + depth);
    output.classList.remove("string_json", "margin_left_30", "max_width_500");
    output.innerHTML = "";
    let ul = document.createElement("ul");
    ul.classList.add("nodes");
    ul.id = "nodes_" + depth;
    output.appendChild(ul);
    let div = document.createElement("div");
    div.classList.add("flex");
    div.id = "nested_div_" + (depth + 1);
    output.appendChild(div);
    config[category].forEach(d => {
      let li = document.createElement("li");
      li.classList.add("node", "node_" + depth);
      li.category = d.id;
      li.id = "node_" + depth + "_" + d.id;
      li.innerHTML = id2label[d.id];
      if (d.api) li.api = d.api;
      ul.appendChild(li);
      if (depth == 0 && !categories.includes(d.id)) return false;
      li.classList.add("active");
      li.addEventListener("click", e => {
	document.querySelectorAll(".node_" + depth).forEach(el => {
	  el.classList.remove("select");
	});
	e.target.classList.add("select");
	if (depth == 0) {
	  showCategoryList(depth + 1, e.target.category, false, false);
	} else {
	  showLocationList(depth + 1, e.target.api, e.target.category, query);
	}
      });
    });
  }


  const showLocationList = async (depth, api_name, selected_category, query) => {
    if (!query) query = document.querySelector("#query").value.replace(/\s/g, "");
    let [, id, location ] = query.match(/^([^:]+):(.+)$/);
    location = location.replace(/(\d)c\d+/g, '$1');
    let category;
    if (id.match(/#/)) [category, id] = id.split(/#/);
    const output = document.querySelector("#nested_div_" + depth);
    output.innerHTML = "";
    let img = document.createElement("img");
    img.classList.add("loading_icon");
    img.src = "./assets/loading.gif";
    output.appendChild(img); 
    options.body = "source=" + id + "&location=" + encodeURIComponent(location);
    await fetch(api + api_name, options).then(r => r.json()).then(json => {
      output.innerHTML = "";
      if (!json[0]) {
	output.innerHTML = "<p class='margin_left_30'>No corresponding coordinates found</p>";
	return false;
      }
      let ul = document.createElement("ul");
      ul.classList.add("locations");
      ul.id = "locations_" + depth;
      output.appendChild(ul);
      let div = document.createElement("div");
      div.classList.add("flex");
      div.id = "nested_div_" + (depth + 1);
      output.appendChild(div);
      json.forEach(d => {
	let li = document.createElement("li");
	li.classList.add("location", "location_" + depth);
	li.location = d.jsonld.id.split(/\//).at(-1);
	li.id = "locatio_" + depth + "_" + li.location.split(/:/)[0];
	li.innerHTML = "<span class='bold'>Converted location:</span><br>";
	li.innerHTML += li.location;
	li.innerHTML += "<br><button class='small_button small_button_active'>JSON-LD</button><button class='small_button small_button_active'>TABLE</button><button class='small_button small_button_inactive' category='" + selected_category + "'>CONVERT</button>";
	li.pre_api = api_name;
	ul.appendChild(li);
	li.querySelectorAll(".small_button")[0].addEventListener("click", e => {
	  showJson(depth + 1, e.target, d.jsonld);
	});
	li.querySelectorAll(".small_button")[1].addEventListener("click", e => {
	  showTable(depth + 1, e.target, d.correspondence);
	});
	if (config[selected_category] && !li.location.match(/complement/)) {
	  let button = li.querySelectorAll(".small_button")[2];
	  button.classList.remove("small_button_inactive");
	  button.classList.add("small_button_active");
	  button.addEventListener("click", e => {
	    showCategoryList(depth + 1, e.target.getAttribute("category"), li.location, false);
	  });
	}
      });
    });
  }
  
  const showJson = (depth, element, json) => {
    element.parentNode.parentNode.querySelectorAll(".location").forEach(element => {
      element.classList.remove("location_select");
    });
    element.parentNode.classList.add("location_select");
    const output = document.querySelector("#nested_div_" + depth);
    output.classList.add("string_json", "margin_left_30", "max_width_500");
    output.innerHTML = "";
    let pre = document.createElement("pre");
    pre.classList.add("jsonld");
    pre.innerHTML =  JSON.stringify(json, null, 2);
    output.appendChild(pre);
  }

  const showTable = (depth, element, json) => {
    element.parentNode.parentNode.querySelectorAll(".location").forEach(element => {
      element.classList.remove("location_select");
    });
    element.parentNode.classList.add("location_select");
    const output = document.querySelector("#nested_div_" + depth);
    output.classList.add("string_json", "margin_left_30", "max_width_500");
    output.innerHTML = "<div class='w500_dummy_for_shadow_root'></div>";
    // tsv in <pre>
   /* let pre = document.createElement("pre");
    let tsv = "Source\tTarget\n"
    for (let d of json.ranges) {
      tsv += d.source_begin + ".." + d.source_end + "\t";
      if (json.complement) tsv += "complement(";
      tsv += d.begin + ".." + d.end;
      if (json.complement) tsv += ")";
      tsv += "\n";
    }
    pre.innerHTML =  tsv;
    output.appendChild(pre); */
    json.ranges.forEach((d, i) => {
      if (!d.joined_complement) json.ranges[i].joined_complement = false;
    });
    const blob = new Blob([JSON.stringify(json.ranges)], { type: 'text/plain' });
    const dataUrl = URL.createObjectURL(blob)
    let stanza = document.createElement("togostanza-pagination-table");
    stanza.setAttribute("data-url", dataUrl);
    stanza.setAttribute("data-type", "json");
    stanza.setAttribute("width", 500);
    stanza.setAttribute("page-size-option", "10,50");
    stanza.setAttribute("columns", '[{"id":"source_begin","label":"Source begin"},{"id":"source_end","label":"Source end"},{"id":"begin","label":"Target begin"},{"id":"end","label":"Target end"},{"id":"joined_complement","label":"Complement","type":"boolean"}]');
    output.appendChild(stanza);
  }
  
  const querySubmit = async () => {
    let query = document.querySelector("#query").value.replace(/\s/g, "");
    if (!query || !query.match(/.+:.*\d/)) return false;
    let categories = [];
    config = await fetch("./assets/config.json").then(r => r.json());
    [query, categories[0]] = checkQueryCategory(query, config);
    const mode = checkMode();
    
    if (mode == "convert") {
      id2label = {};
      config.categories.forEach(d => {
	id2label[d.id] = d.label;
      });
      if (!categories[0]) categories = getQueryCategories(query);
      showCategoryList(0, "all", false, categories);
 
    } else if (mode == "jsonld") {
      let prefix = document.querySelector("#prefix").value || "undefined";
      options.body = "location_id=" + encodeURIComponent(query) + "&reference_prefix=" + prefix;
      await fetch(api + "togocoord_location_converter", options).then(r => r.json()).then(j => {
	document.querySelector("#output_jsonld").innerHTML = JSON.stringify(j, null, 2);
      });
    }
  }

  const queryReset = () => {
    document.querySelector("#query").value = "";
    document.querySelector("#prefix").value = "";
    document.querySelector("#nested_div_0").innerHTML = "";
    document.querySelector("#output_jsonld").innerHTML = "";
  }

  //// Event listener
  // Query
  document.querySelector("#submit").addEventListener("click", querySubmit);
  document.querySelector("#reset").addEventListener("click", queryReset);
  // Example
  document.querySelectorAll("button.example").forEach(element => {
    element.addEventListener("click", event => {
      queryReset();
      document.querySelector("#query").value = event.target.getAttribute("location");
      querySubmit();
    });
  });
  // Tab
  document.querySelectorAll(".tab").forEach(element => {
    element.addEventListener("click", event => {
      document.querySelectorAll(".tab").forEach(element => {
	element.classList.remove("tab_selected");
	element.classList.add("tab_unselected");
      });
      event.target.classList.remove("tab_unselected");
      event.target.classList.add("tab_selected");
      const id = "#result_" + event.target.id;
      document.querySelectorAll(".result").forEach(element => {
	element.classList.remove("show");
	element.classList.add("hidden");
      });
      document.querySelector(id).classList.remove("hidden");
      document.querySelector(id).classList.add("show");
    });
  });

  //// Auto scroll
  const observer = new MutationObserver((mutationsList) => {
    for (let mutation of mutationsList) {
      if (mutation.type === 'childList') {
	requestAnimationFrame(() => {
          requestAnimationFrame(() => {
	    let element = document.querySelector("#result");
	    console.log([element.scrollLeft, element.scrollWidth]);
	    element.scrollLeft = element.scrollWidth;
          });
	});
      }
    }
  });
  observer.observe(document.body, { childList: true, subtree: true });
  
};

TogoCoord();
