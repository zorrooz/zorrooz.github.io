const __vite__mapDeps=(i,m=__vite__mapDeps,d=(m.f||(m.f=["assets/SearchModal-CPATOImV.js","assets/SearchModal-DvrNkWgZ.css"])))=>i.map(i=>d[i]);
(function(){const n=document.createElement("link").relList;if(n&&n.supports&&n.supports("modulepreload"))return;for(const t of document.querySelectorAll('link[rel="modulepreload"]'))e(t);new MutationObserver(t=>{for(const l of t)if(l.type==="childList")for(const o of l.addedNodes)o.tagName==="LINK"&&o.rel==="modulepreload"&&e(o)}).observe(document,{childList:!0,subtree:!0});function a(t){const l={};return t.integrity&&(l.integrity=t.integrity),t.referrerPolicy&&(l.referrerPolicy=t.referrerPolicy),t.crossOrigin==="use-credentials"?l.credentials="include":t.crossOrigin==="anonymous"?l.credentials="omit":l.credentials="same-origin",l}function e(t){if(t.ep)return;t.ep=!0;const l=a(t);fetch(t.href,l)}})();const Ji="modulepreload",Zi=function(s){return"/"+s},xo={},ys=function(n,a,e){let t=Promise.resolve();if(a&&a.length>0){let c=function(u){return Promise.all(u.map(r=>Promise.resolve(r).then(i=>({status:"fulfilled",value:i}),i=>({status:"rejected",reason:i}))))};document.getElementsByTagName("link");const o=document.querySelector("meta[property=csp-nonce]"),p=o?.nonce||o?.getAttribute("nonce");t=c(a.map(u=>{if(u=Zi(u),u in xo)return;xo[u]=!0;const r=u.endsWith(".css"),i=r?'[rel="stylesheet"]':"";if(document.querySelector(`link[href="${u}"]${i}`))return;const d=document.createElement("link");if(d.rel=r?"stylesheet":Ji,r||(d.as="script"),d.crossOrigin="",d.href=u,p&&d.setAttribute("nonce",p),document.head.appendChild(d),r)return new Promise((h,y)=>{d.addEventListener("load",h),d.addEventListener("error",()=>y(new Error(`Unable to preload CSS for ${u}`)))})}))}function l(o){const p=new Event("vite:preloadError",{cancelable:!0});if(p.payload=o,window.dispatchEvent(p),!p.defaultPrevented)throw o}return t.then(o=>{for(const p of o||[])p.status==="rejected"&&l(p.reason);return n().catch(l)})};function pl(s,n={},a){for(const e in s){const t=s[e],l=a?`${a}:${e}`:e;typeof t=="object"&&t!==null?pl(t,n,l):typeof t=="function"&&(n[l]=t)}return n}const yc=(()=>{if(console.createTask)return console.createTask;const s={run:n=>n()};return()=>s})();function vc(s,n,a,e){for(let t=a;t<s.length;t+=1)try{const l=e?e.run(()=>s[t](...n)):s[t](...n);if(l&&typeof l.then=="function")return Promise.resolve(l).then(()=>vc(s,n,t+1,e))}catch(l){return Promise.reject(l)}}function su(s,n,a){if(s.length>0)return vc(s,n,0,yc(a))}function nu(s,n,a){if(s.length>0){const e=yc(a);return Promise.all(s.map(t=>e.run(()=>t(...n))))}}function $t(s,n){for(const a of[...s])a(n)}var au=class{_hooks;_before;_after;_deprecatedHooks;_deprecatedMessages;constructor(){this._hooks={},this._before=void 0,this._after=void 0,this._deprecatedMessages=void 0,this._deprecatedHooks={},this.hook=this.hook.bind(this),this.callHook=this.callHook.bind(this),this.callHookWith=this.callHookWith.bind(this)}hook(s,n,a={}){if(!s||typeof n!="function")return()=>{};const e=s;let t;for(;this._deprecatedHooks[s];)t=this._deprecatedHooks[s],s=t.to;if(t&&!a.allowDeprecated){let l=t.message;l||(l=`${e} hook has been deprecated`+(t.to?`, please use ${t.to}`:"")),this._deprecatedMessages||(this._deprecatedMessages=new Set),this._deprecatedMessages.has(l)||(console.warn(l),this._deprecatedMessages.add(l))}if(!n.name)try{Object.defineProperty(n,"name",{get:()=>"_"+s.replace(/\W+/g,"_")+"_hook_cb",configurable:!0})}catch{}return this._hooks[s]=this._hooks[s]||[],this._hooks[s].push(n),()=>{n&&(this.removeHook(s,n),n=void 0)}}hookOnce(s,n){let a,e=(...t)=>(typeof a=="function"&&a(),a=void 0,e=void 0,n(...t));return a=this.hook(s,e),a}removeHook(s,n){const a=this._hooks[s];if(a){const e=a.indexOf(n);e!==-1&&a.splice(e,1),a.length===0&&(this._hooks[s]=void 0)}}clearHook(s){this._hooks[s]=void 0}deprecateHook(s,n){this._deprecatedHooks[s]=typeof n=="string"?{to:n}:n;const a=this._hooks[s]||[];this._hooks[s]=void 0;for(const e of a)this.hook(s,e)}deprecateHooks(s){for(const n in s)this.deprecateHook(n,s[n])}addHooks(s){const n=pl(s),a=Object.keys(n).map(e=>this.hook(e,n[e]));return()=>{for(const e of a)e();a.length=0}}removeHooks(s){const n=pl(s);for(const a in n)this.removeHook(a,n[a])}removeAllHooks(){this._hooks={}}callHook(s,...n){return this.callHookWith(su,s,n)}callHookParallel(s,...n){return this.callHookWith(nu,s,n)}callHookWith(s,n,a){const e=this._before||this._after?{name:n,args:a,context:{}}:void 0;this._before&&$t(this._before,e);const t=s(this._hooks[n]?[...this._hooks[n]]:[],a,n);return t instanceof Promise?t.finally(()=>{this._after&&e&&$t(this._after,e)}):(this._after&&e&&$t(this._after,e),t)}beforeEach(s){return this._before=this._before||[],this._before.push(s),()=>{if(this._before!==void 0){const n=this._before.indexOf(s);n!==-1&&this._before.splice(n,1)}}}afterEach(s){return this._after=this._after||[],this._after.push(s),()=>{if(this._after!==void 0){const n=this._after.indexOf(s);n!==-1&&this._after.splice(n,1)}}}};function eu(){return new au}const tu=new Set(["link","style","script","noscript"]),lu=new Set(["title","titleTemplate","script","style","noscript"]),cl=new Set(["base","meta","link","style","script","noscript"]),ou=new Set(["title","base","htmlAttrs","bodyAttrs","meta","link","style","script","noscript"]),pu=new Set(["base","title","titleTemplate","bodyAttrs","htmlAttrs","templateParams"]),cu=new Set(["key","tagPosition","tagPriority","tagDuplicateStrategy","innerHTML","textContent","processTemplateParams"]),ru=new Set(["templateParams","htmlAttrs","bodyAttrs"]),iu=new Set(["theme-color","google-site-verification","og","article","book","profile","twitter","author"]),uu=["name","property","http-equiv"],hu=new Set(["viewport","description","keywords","robots"]);function wc(s){const n=s.split(":");return n.length?iu.has(n[1]):!1}function rl(s){const{props:n,tag:a}=s;if(pu.has(a))return a;if(a==="link"&&n.rel==="canonical")return"canonical";if(a==="link"&&n.rel==="alternate"){if(n.hreflang)return`alternate:${n.hreflang}`;if(n.type)return`alternate:${n.type}:${n.href||""}`}if(n.charset)return"charset";if(s.tag==="meta"){for(const e of uu)if(n[e]!==void 0){const t=n[e],l=t&&typeof t=="string"&&t.includes(":"),o=t&&hu.has(t),c=!(l||o)&&s.key?`:key:${s.key}`:"";return`${a}:${t}${c}`}}if(s.key)return`${a}:key:${s.key}`;if(n.id)return`${a}:id:${n.id}`;if(a==="link"&&n.rel==="alternate")return`alternate:${n.href||""}`;if(lu.has(a)){const e=s.textContent||s.innerHTML;if(e)return`${a}:content:${e}`}}function Cc(s){const n=s._h||s._d;if(n)return n;const a=s.textContent||s.innerHTML;return a||`${s.tag}:${Object.entries(s.props).map(([e,t])=>`${e}:${String(t)}`).join(",")}`}function tt(s,n,a){typeof s==="function"&&(!a||a!=="titleTemplate"&&!(a[0]==="o"&&a[1]==="n"))&&(s=s());const t=n?n(a,s):s;if(Array.isArray(t))return t.map(l=>tt(l,n));if(t?.constructor===Object){const l={};for(const o of Object.keys(t))l[o]=tt(t[o],n,o);return l}return t}function du(s,n){const a=s==="style"?new Map:new Set;function e(t){if(t==null||t===void 0)return;const l=String(t).trim();if(l)if(s==="style"){const[o,...p]=l.split(":").map(c=>c?c.trim():"");o&&p.length&&a.set(o,p.join(":"))}else l.split(" ").filter(Boolean).forEach(o=>a.add(o))}return typeof n=="string"?s==="style"?n.split(";").forEach(e):e(n):Array.isArray(n)?n.forEach(t=>e(t)):n&&typeof n=="object"&&Object.entries(n).forEach(([t,l])=>{l&&l!=="false"&&(s==="style"?a.set(String(t).trim(),String(l)):e(t))}),a}function kc(s,n){if(s.props=s.props||{},!n)return s;if(s.tag==="templateParams")return s.props=n,s;const a=cl.has(s.tag)||s.tag==="htmlAttrs"||s.tag==="bodyAttrs";return Object.entries(n).forEach(([e,t])=>{if(e==="__proto__"||e==="constructor"||e==="prototype")return;if(t===null){s.props[e]=null;return}if(e==="class"||e==="style"){s.props[e]=du(e,t);return}if(cu.has(e)){if((e==="textContent"||e==="innerHTML")&&typeof t=="object"){let u=n.type;if(n.type||(u="application/json"),!u?.endsWith("json")&&u!=="speculationrules")return;n.type=u,s.props.type=u,s[e]=JSON.stringify(t)}else s[e]=t;return}const l=e.startsWith("data-"),o=a&&!l?e.toLowerCase():e,p=String(t),c=s.tag==="meta"&&o==="content";p==="true"||p===""?s.props[o]=l||c?p:!0:!t&&l&&p==="false"?s.props[o]="false":t!==void 0&&(s.props[o]=t)}),s}function mu(s,n){const a=typeof n=="object"&&typeof n!="function"?n:{[s==="script"||s==="noscript"||s==="style"?"innerHTML":"textContent"]:n},e=kc({tag:s,props:{}},a);return e.key&&tu.has(e.tag)&&(e.props["data-hid"]=e._h=e.key),e.tag==="script"&&typeof e.innerHTML=="object"&&(e.innerHTML=JSON.stringify(e.innerHTML),e.props.type=e.props.type||"application/json"),Array.isArray(e.props.content)?e.props.content.map(t=>({...e,props:{...e.props,content:t}})):e}function ju(s,n){if(!s)return[];typeof s=="function"&&(s=s());const a=(t,l)=>{for(let o=0;o<n.length;o++)l=n[o](t,l);return l};s=a(void 0,s);const e=[];return s=tt(s,a),Object.entries(s||{}).forEach(([t,l])=>{if(l!==void 0)for(const o of Array.isArray(l)?l:[l])e.push(mu(t,o))}),e.flat()}const Po=(s,n)=>s._w===n._w?s._p-n._p:s._w-n._w,Eo={base:-10,title:10},fu={critical:-8,high:-1,low:2},So={meta:{"content-security-policy":-30,charset:-20,viewport:-15},link:{preconnect:20,stylesheet:60,preload:70,modulepreload:70,prefetch:90,"dns-prefetch":90,prerender:90},script:{async:30,defer:80,sync:50},style:{imported:40,sync:60}},gu=/@import/,Ya=s=>s===""||s===!0;function bu(s,n){if(typeof n.tagPriority=="number")return n.tagPriority;let a=100;const e=fu[n.tagPriority]||0,t=s.resolvedOptions.disableCapoSorting?{link:{},script:{},style:{}}:So;if(n.tag in Eo)a=Eo[n.tag];else if(n.tag==="meta"){const l=n.props["http-equiv"]==="content-security-policy"?"content-security-policy":n.props.charset?"charset":n.props.name==="viewport"?"viewport":null;l&&(a=So.meta[l])}else if(n.tag==="link"&&n.props.rel)a=t.link[n.props.rel];else if(n.tag==="script"){const l=String(n.props.type);Ya(n.props.async)?a=t.script.async:n.props.src&&!Ya(n.props.defer)&&!Ya(n.props.async)&&l!=="module"&&!l.endsWith("json")||n.innerHTML&&!l.endsWith("json")?a=t.script.sync:(Ya(n.props.defer)&&n.props.src&&!Ya(n.props.async)||l==="module")&&(a=t.script.defer)}else n.tag==="style"&&(a=n.innerHTML&&gu.test(n.innerHTML)?t.style.imported:t.style.sync);return(a||100)+e}function To(s,n){const a=typeof n=="function"?n(s):n,e=a.key||String(s.plugins.size+1);s.plugins.get(e)||(s.plugins.set(e,a),s.hooks.addHooks(a.hooks||{}))}function _u(s={}){const n=eu();n.addHooks(s.hooks||{});const a=!s.document,e=new Map,t=new Map,l=new Set,o={_entryCount:1,plugins:t,dirty:!1,resolvedOptions:s,hooks:n,ssr:a,entries:e,headEntries(){return[...e.values()]},use:p=>To(o,p),push(p,c){const u={...c||{}};delete u.head;const r=u._index??o._entryCount++,i={_i:r,input:p,options:u},d={_poll(h=!1){o.dirty=!0,!h&&l.add(r),n.callHook("entries:updated",o)},dispose(){e.delete(r)&&o.invalidate()},patch(h){(!u.mode||u.mode==="server"&&a||u.mode==="client"&&!a)&&(i.input=h,e.set(r,i),d._poll())}};return d.patch(p),d},async resolveTags(){const p={tagMap:new Map,tags:[],entries:[...o.entries.values()]};for(await n.callHook("entries:resolve",p);l.size;){const d=l.values().next().value;l.delete(d);const h=e.get(d);if(h){const y={tags:ju(h.input,s.propResolvers||[]).map(j=>Object.assign(j,h.options)),entry:h};await n.callHook("entries:normalize",y),h._tags=y.tags.map((j,T)=>(j._w=bu(o,j),j._p=(h._i<<10)+T,j._d=rl(j),j._d||(j._h=Cc(j)),j))}}let c=!1;p.entries.flatMap(d=>(d._tags||[]).map(h=>({...h,props:{...h.props}}))).sort(Po).reduce((d,h)=>{const y=h._d||h._h;if(!d.has(y))return d.set(y,h);const j=d.get(y);if((h?.tagDuplicateStrategy||(ru.has(h.tag)?"merge":null)||(h.key&&h.key===j.key?"merge":null))==="merge"){const g={...j.props};Object.entries(h.props).forEach(([m,w])=>g[m]=m==="style"?new Map([...j.props.style||new Map,...w]):m==="class"?new Set([...j.props.class||new Set,...w]):w),d.set(y,{...h,props:g})}else h._p>>10===j._p>>10&&h.tag==="meta"&&wc(y)?(d.set(y,Object.assign([...Array.isArray(j)?j:[j],h],h)),c=!0):(h._w===j._w?h._p>j._p:h?._w<j?._w)&&d.set(y,h);return d},p.tagMap);const u=p.tagMap.get("title"),r=p.tagMap.get("titleTemplate");if(o._title=u?.textContent,r){const d=r?.textContent;if(o._titleTemplate=d,d){let h=typeof d=="function"?d(u?.textContent):d;typeof h=="string"&&!o.plugins.has("template-params")&&(h=h.replace("%s",u?.textContent||"")),u?h===null?p.tagMap.delete("title"):p.tagMap.set("title",{...u,textContent:h}):(r.tag="title",r.textContent=h)}}p.tags=Array.from(p.tagMap.values()),c&&(p.tags=p.tags.flat().sort(Po)),await n.callHook("tags:beforeResolve",p),await n.callHook("tags:resolve",p),await n.callHook("tags:afterResolve",p);const i=[];for(const d of p.tags){const{innerHTML:h,tag:y,props:j}=d;if(ou.has(y)&&!(Object.keys(j).length===0&&!d.innerHTML&&!d.textContent)&&!(y==="meta"&&!j.content&&!j["http-equiv"]&&!j.charset)){if(y==="script"&&h){if(String(j.type).endsWith("json")){const T=typeof h=="string"?h:JSON.stringify(h);d.innerHTML=T.replace(/</g,"\\u003C")}else typeof h=="string"&&(d.innerHTML=h.replace(new RegExp(`</${y}`,"g"),`<\\/${y}`));d._d=rl(d)}i.push(d)}}return i},invalidate(){for(const p of e.values())l.add(p._i);o.dirty=!0,n.callHook("entries:updated",o)}};return(s?.plugins||[]).forEach(p=>To(o,p)),o.hooks.callHook("init",o),s.init?.forEach(p=>p&&o.push(p)),o}async function xc(s,n={}){const a=n.document||s.resolvedOptions.document;if(!a||!s.dirty)return;const e={shouldRender:!0,tags:[]};if(await s.hooks.callHook("dom:beforeRender",e),!!e.shouldRender)return s._domUpdatePromise||(s._domUpdatePromise=new Promise(async t=>{const l=new Map,o=new Promise(h=>{s.resolveTags().then(y=>{h(y.map(j=>{const T=l.get(j._d)||0,g={tag:j,id:(T?`${j._d}:${T}`:j._d)||j._h,shouldRender:!0};return j._d&&wc(j._d)&&l.set(j._d,T+1),g}))})});let p=s._dom;if(!p){p={title:a.title,elMap:new Map().set("htmlAttrs",a.documentElement).set("bodyAttrs",a.body)};for(const h of["body","head"]){const y=a[h]?.children;for(const j of y){const T=j.tagName.toLowerCase();if(!cl.has(T))continue;const g=kc({tag:T,props:{}},{innerHTML:j.innerHTML,...j.getAttributeNames().reduce((m,w)=>(m[w]=j.getAttribute(w),m),{})||{}});if(g.key=j.getAttribute("data-hid")||void 0,g._d=rl(g)||Cc(g),p.elMap.has(g._d)){let m=1,w=g._d;for(;p.elMap.has(w);)w=`${g._d}:${m++}`;p.elMap.set(w,j)}else p.elMap.set(g._d,j)}}}p.pendingSideEffects={...p.sideEffects},p.sideEffects={};function c(h,y,j){const T=`${h}:${y}`;p.sideEffects[T]=j,delete p.pendingSideEffects[T]}function u({id:h,$el:y,tag:j}){const T=j.tag.endsWith("Attrs");p.elMap.set(h,y),T||(j.textContent&&j.textContent!==y.textContent&&(y.textContent=j.textContent),j.innerHTML&&j.innerHTML!==y.innerHTML&&(y.innerHTML=j.innerHTML),c(h,"el",()=>{y?.remove(),p.elMap.delete(h)}));for(const g in j.props){if(!Object.prototype.hasOwnProperty.call(j.props,g))continue;const m=j.props[g];if(g.startsWith("on")&&typeof m=="function"){const f=y?.dataset;if(f&&f[`${g}fired`]){const x=g.slice(0,-5);m.call(y,new Event(x.substring(2)))}y.getAttribute(`data-${g}`)!==""&&((j.tag==="bodyAttrs"?a.defaultView:y).addEventListener(g.substring(2),m.bind(y)),y.setAttribute(`data-${g}`,""));continue}const w=`attr:${g}`;if(g==="class"){if(!m)continue;for(const f of m)T&&c(h,`${w}:${f}`,()=>y.classList.remove(f)),!y.classList.contains(f)&&y.classList.add(f)}else if(g==="style"){if(!m)continue;for(const[f,x]of m)c(h,`${w}:${f}`,()=>{y.style.removeProperty(f)}),y.style.setProperty(f,x)}else m!==!1&&m!==null&&(y.getAttribute(g)!==m&&y.setAttribute(g,m===!0?"":String(m)),T&&c(h,w,()=>y.removeAttribute(g)))}}const r=[],i={bodyClose:void 0,bodyOpen:void 0,head:void 0},d=await o;for(const h of d){const{tag:y,shouldRender:j,id:T}=h;if(j){if(y.tag==="title"){a.title=y.textContent,c("title","",()=>a.title=p.title);continue}h.$el=h.$el||p.elMap.get(T),h.$el?u(h):cl.has(y.tag)&&r.push(h)}}for(const h of r){const y=h.tag.tagPosition||"head";h.$el=a.createElement(h.tag.tag),u(h),i[y]=i[y]||a.createDocumentFragment(),i[y].appendChild(h.$el)}for(const h of d)await s.hooks.callHook("dom:renderTag",h,a,c);i.head&&a.head.appendChild(i.head),i.bodyOpen&&a.body.insertBefore(i.bodyOpen,a.body.firstChild),i.bodyClose&&a.body.appendChild(i.bodyClose);for(const h in p.pendingSideEffects)p.pendingSideEffects[h]();s._dom=p,await s.hooks.callHook("dom:rendered",{renders:d}),t()}).finally(()=>{s._domUpdatePromise=void 0,s.dirty=!1})),s._domUpdatePromise}function yu(s={}){const n=s.domOptions?.render||xc;s.document=s.document||(typeof window<"u"?document:void 0);const a=s.document?.head.querySelector('script[id="unhead:payload"]')?.innerHTML||!1;return _u({...s,plugins:[...s.plugins||[],{key:"client",hooks:{"entries:updated":n}}],init:[a?JSON.parse(a):!1,...s.init||[]]})}function vu(s,n){let a=0;return()=>{const e=++a;n(()=>{a===e&&s()})}}/**
* @vue/shared v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/function Bl(s){const n=Object.create(null);for(const a of s.split(","))n[a]=1;return a=>a in n}const Ss={},Ba=[],$n=()=>{},Pc=()=>!1,yt=s=>s.charCodeAt(0)===111&&s.charCodeAt(1)===110&&(s.charCodeAt(2)>122||s.charCodeAt(2)<97),ql=s=>s.startsWith("onUpdate:"),Ks=Object.assign,zl=(s,n)=>{const a=s.indexOf(n);a>-1&&s.splice(a,1)},wu=Object.prototype.hasOwnProperty,Es=(s,n)=>wu.call(s,n),ds=Array.isArray,qa=s=>vt(s)==="[object Map]",Ec=s=>vt(s)==="[object Set]",js=s=>typeof s=="function",$s=s=>typeof s=="string",fa=s=>typeof s=="symbol",Ds=s=>s!==null&&typeof s=="object",Sc=s=>(Ds(s)||js(s))&&js(s.then)&&js(s.catch),Tc=Object.prototype.toString,vt=s=>Tc.call(s),Cu=s=>vt(s).slice(8,-1),Ac=s=>vt(s)==="[object Object]",Vl=s=>$s(s)&&s!=="NaN"&&s[0]!=="-"&&""+parseInt(s,10)===s,te=Bl(",key,ref,ref_for,ref_key,onVnodeBeforeMount,onVnodeMounted,onVnodeBeforeUpdate,onVnodeUpdated,onVnodeBeforeUnmount,onVnodeUnmounted"),wt=s=>{const n=Object.create(null);return(a=>n[a]||(n[a]=s(a)))},ku=/-\w/g,Tn=wt(s=>s.replace(ku,n=>n.slice(1).toUpperCase())),xu=/\B([A-Z])/g,ga=wt(s=>s.replace(xu,"-$1").toLowerCase()),Ct=wt(s=>s.charAt(0).toUpperCase()+s.slice(1)),Bt=wt(s=>s?`on${Ct(s)}`:""),ha=(s,n)=>!Object.is(s,n),Ye=(s,...n)=>{for(let a=0;a<s.length;a++)s[a](...n)},Rc=(s,n,a,e=!1)=>{Object.defineProperty(s,n,{configurable:!0,enumerable:!1,writable:e,value:a})},Ul=s=>{const n=parseFloat(s);return isNaN(n)?s:n},Pu=s=>{const n=$s(s)?Number(s):NaN;return isNaN(n)?s:n};let Ao;const kt=()=>Ao||(Ao=typeof globalThis<"u"?globalThis:typeof self<"u"?self:typeof window<"u"?window:typeof global<"u"?global:{});function ba(s){if(ds(s)){const n={};for(let a=0;a<s.length;a++){const e=s[a],t=$s(e)?Au(e):ba(e);if(t)for(const l in t)n[l]=t[l]}return n}else if($s(s)||Ds(s))return s}const Eu=/;(?![^(]*\))/g,Su=/:([^]+)/,Tu=/\/\*[^]*?\*\//g;function Au(s){const n={};return s.replace(Tu,"").split(Eu).forEach(a=>{if(a){const e=a.split(Su);e.length>1&&(n[e[0].trim()]=e[1].trim())}}),n}function Us(s){let n="";if($s(s))n=s;else if(ds(s))for(let a=0;a<s.length;a++){const e=Us(s[a]);e&&(n+=e+" ")}else if(Ds(s))for(const a in s)s[a]&&(n+=a+" ");return n.trim()}const Ru="itemscope,allowfullscreen,formnovalidate,ismap,nomodule,novalidate,readonly",Lu=Bl(Ru);function Lc(s){return!!s||s===""}const Mc=s=>!!(s&&s.__v_isRef===!0),X=s=>$s(s)?s:s==null?"":ds(s)||Ds(s)&&(s.toString===Tc||!js(s.toString))?Mc(s)?X(s.value):JSON.stringify(s,Dc,2):String(s),Dc=(s,n)=>Mc(n)?Dc(s,n.value):qa(n)?{[`Map(${n.size})`]:[...n.entries()].reduce((a,[e,t],l)=>(a[qt(e,l)+" =>"]=t,a),{})}:Ec(n)?{[`Set(${n.size})`]:[...n.values()].map(a=>qt(a))}:fa(n)?qt(n):Ds(n)&&!ds(n)&&!Ac(n)?String(n):n,qt=(s,n="")=>{var a;return fa(s)?`Symbol(${(a=s.description)!=null?a:n})`:s};/**
* @vue/reactivity v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/let en;class Oc{constructor(n=!1){this.detached=n,this._active=!0,this._on=0,this.effects=[],this.cleanups=[],this._isPaused=!1,this.parent=en,!n&&en&&(this.index=(en.scopes||(en.scopes=[])).push(this)-1)}get active(){return this._active}pause(){if(this._active){this._isPaused=!0;let n,a;if(this.scopes)for(n=0,a=this.scopes.length;n<a;n++)this.scopes[n].pause();for(n=0,a=this.effects.length;n<a;n++)this.effects[n].pause()}}resume(){if(this._active&&this._isPaused){this._isPaused=!1;let n,a;if(this.scopes)for(n=0,a=this.scopes.length;n<a;n++)this.scopes[n].resume();for(n=0,a=this.effects.length;n<a;n++)this.effects[n].resume()}}run(n){if(this._active){const a=en;try{return en=this,n()}finally{en=a}}}on(){++this._on===1&&(this.prevScope=en,en=this)}off(){this._on>0&&--this._on===0&&(en=this.prevScope,this.prevScope=void 0)}stop(n){if(this._active){this._active=!1;let a,e;for(a=0,e=this.effects.length;a<e;a++)this.effects[a].stop();for(this.effects.length=0,a=0,e=this.cleanups.length;a<e;a++)this.cleanups[a]();if(this.cleanups.length=0,this.scopes){for(a=0,e=this.scopes.length;a<e;a++)this.scopes[a].stop(!0);this.scopes.length=0}if(!this.detached&&this.parent&&!n){const t=this.parent.scopes.pop();t&&t!==this&&(this.parent.scopes[this.index]=t,t.index=this.index)}this.parent=void 0}}}function Hl(s){return new Oc(s)}function Ic(){return en}function Mu(s,n=!1){en&&en.cleanups.push(s)}let Ls;const zt=new WeakSet;class Nc{constructor(n){this.fn=n,this.deps=void 0,this.depsTail=void 0,this.flags=5,this.next=void 0,this.cleanup=void 0,this.scheduler=void 0,en&&en.active&&en.effects.push(this)}pause(){this.flags|=64}resume(){this.flags&64&&(this.flags&=-65,zt.has(this)&&(zt.delete(this),this.trigger()))}notify(){this.flags&2&&!(this.flags&32)||this.flags&8||$c(this)}run(){if(!(this.flags&1))return this.fn();this.flags|=2,Ro(this),Bc(this);const n=Ls,a=Rn;Ls=this,Rn=!0;try{return this.fn()}finally{qc(this),Ls=n,Rn=a,this.flags&=-3}}stop(){if(this.flags&1){for(let n=this.deps;n;n=n.nextDep)Kl(n);this.deps=this.depsTail=void 0,Ro(this),this.onStop&&this.onStop(),this.flags&=-2}}trigger(){this.flags&64?zt.add(this):this.scheduler?this.scheduler():this.runIfDirty()}runIfDirty(){il(this)&&this.run()}get dirty(){return il(this)}}let Fc=0,le,oe;function $c(s,n=!1){if(s.flags|=8,n){s.next=oe,oe=s;return}s.next=le,le=s}function Gl(){Fc++}function Wl(){if(--Fc>0)return;if(oe){let n=oe;for(oe=void 0;n;){const a=n.next;n.next=void 0,n.flags&=-9,n=a}}let s;for(;le;){let n=le;for(le=void 0;n;){const a=n.next;if(n.next=void 0,n.flags&=-9,n.flags&1)try{n.trigger()}catch(e){s||(s=e)}n=a}}if(s)throw s}function Bc(s){for(let n=s.deps;n;n=n.nextDep)n.version=-1,n.prevActiveLink=n.dep.activeLink,n.dep.activeLink=n}function qc(s){let n,a=s.depsTail,e=a;for(;e;){const t=e.prevDep;e.version===-1?(e===a&&(a=t),Kl(e),Du(e)):n=e,e.dep.activeLink=e.prevActiveLink,e.prevActiveLink=void 0,e=t}s.deps=n,s.depsTail=a}function il(s){for(let n=s.deps;n;n=n.nextDep)if(n.dep.version!==n.version||n.dep.computed&&(zc(n.dep.computed)||n.dep.version!==n.version))return!0;return!!s._dirty}function zc(s){if(s.flags&4&&!(s.flags&16)||(s.flags&=-17,s.globalVersion===fe)||(s.globalVersion=fe,!s.isSSR&&s.flags&128&&(!s.deps&&!s._dirty||!il(s))))return;s.flags|=2;const n=s.dep,a=Ls,e=Rn;Ls=s,Rn=!0;try{Bc(s);const t=s.fn(s._value);(n.version===0||ha(t,s._value))&&(s.flags|=128,s._value=t,n.version++)}catch(t){throw n.version++,t}finally{Ls=a,Rn=e,qc(s),s.flags&=-3}}function Kl(s,n=!1){const{dep:a,prevSub:e,nextSub:t}=s;if(e&&(e.nextSub=t,s.prevSub=void 0),t&&(t.prevSub=e,s.nextSub=void 0),a.subs===s&&(a.subs=e,!e&&a.computed)){a.computed.flags&=-5;for(let l=a.computed.deps;l;l=l.nextDep)Kl(l,!0)}!n&&!--a.sc&&a.map&&a.map.delete(a.key)}function Du(s){const{prevDep:n,nextDep:a}=s;n&&(n.nextDep=a,s.prevDep=void 0),a&&(a.prevDep=n,s.nextDep=void 0)}let Rn=!0;const Vc=[];function Zn(){Vc.push(Rn),Rn=!1}function sa(){const s=Vc.pop();Rn=s===void 0?!0:s}function Ro(s){const{cleanup:n}=s;if(s.cleanup=void 0,n){const a=Ls;Ls=void 0;try{n()}finally{Ls=a}}}let fe=0;class Ou{constructor(n,a){this.sub=n,this.dep=a,this.version=a.version,this.nextDep=this.prevDep=this.nextSub=this.prevSub=this.prevActiveLink=void 0}}class Xl{constructor(n){this.computed=n,this.version=0,this.activeLink=void 0,this.subs=void 0,this.map=void 0,this.key=void 0,this.sc=0,this.__v_skip=!0}track(n){if(!Ls||!Rn||Ls===this.computed)return;let a=this.activeLink;if(a===void 0||a.sub!==Ls)a=this.activeLink=new Ou(Ls,this),Ls.deps?(a.prevDep=Ls.depsTail,Ls.depsTail.nextDep=a,Ls.depsTail=a):Ls.deps=Ls.depsTail=a,Uc(a);else if(a.version===-1&&(a.version=this.version,a.nextDep)){const e=a.nextDep;e.prevDep=a.prevDep,a.prevDep&&(a.prevDep.nextDep=e),a.prevDep=Ls.depsTail,a.nextDep=void 0,Ls.depsTail.nextDep=a,Ls.depsTail=a,Ls.deps===a&&(Ls.deps=e)}return a}trigger(n){this.version++,fe++,this.notify(n)}notify(n){Gl();try{for(let a=this.subs;a;a=a.prevSub)a.sub.notify()&&a.sub.dep.notify()}finally{Wl()}}}function Uc(s){if(s.dep.sc++,s.sub.flags&4){const n=s.dep.computed;if(n&&!s.dep.subs){n.flags|=20;for(let e=n.deps;e;e=e.nextDep)Uc(e)}const a=s.dep.subs;a!==s&&(s.prevSub=a,a&&(a.nextSub=s)),s.dep.subs=s}}const lt=new WeakMap,Ta=Symbol(""),ul=Symbol(""),ge=Symbol("");function ln(s,n,a){if(Rn&&Ls){let e=lt.get(s);e||lt.set(s,e=new Map);let t=e.get(a);t||(e.set(a,t=new Xl),t.map=e,t.key=a),t.track()}}function Xn(s,n,a,e,t,l){const o=lt.get(s);if(!o){fe++;return}const p=c=>{c&&c.trigger()};if(Gl(),n==="clear")o.forEach(p);else{const c=ds(s),u=c&&Vl(a);if(c&&a==="length"){const r=Number(e);o.forEach((i,d)=>{(d==="length"||d===ge||!fa(d)&&d>=r)&&p(i)})}else switch((a!==void 0||o.has(void 0))&&p(o.get(a)),u&&p(o.get(ge)),n){case"add":c?u&&p(o.get("length")):(p(o.get(Ta)),qa(s)&&p(o.get(ul)));break;case"delete":c||(p(o.get(Ta)),qa(s)&&p(o.get(ul)));break;case"set":qa(s)&&p(o.get(Ta));break}}Wl()}function Iu(s,n){const a=lt.get(s);return a&&a.get(n)}function Ma(s){const n=xs(s);return n===s?n:(ln(n,"iterate",ge),Pn(s)?n:n.map(Qs))}function xt(s){return ln(s=xs(s),"iterate",ge),s}const Nu={__proto__:null,[Symbol.iterator](){return Vt(this,Symbol.iterator,Qs)},concat(...s){return Ma(this).concat(...s.map(n=>ds(n)?Ma(n):n))},entries(){return Vt(this,"entries",s=>(s[1]=Qs(s[1]),s))},every(s,n){return Vn(this,"every",s,n,void 0,arguments)},filter(s,n){return Vn(this,"filter",s,n,a=>a.map(Qs),arguments)},find(s,n){return Vn(this,"find",s,n,Qs,arguments)},findIndex(s,n){return Vn(this,"findIndex",s,n,void 0,arguments)},findLast(s,n){return Vn(this,"findLast",s,n,Qs,arguments)},findLastIndex(s,n){return Vn(this,"findLastIndex",s,n,void 0,arguments)},forEach(s,n){return Vn(this,"forEach",s,n,void 0,arguments)},includes(...s){return Ut(this,"includes",s)},indexOf(...s){return Ut(this,"indexOf",s)},join(s){return Ma(this).join(s)},lastIndexOf(...s){return Ut(this,"lastIndexOf",s)},map(s,n){return Vn(this,"map",s,n,void 0,arguments)},pop(){return Qa(this,"pop")},push(...s){return Qa(this,"push",s)},reduce(s,...n){return Lo(this,"reduce",s,n)},reduceRight(s,...n){return Lo(this,"reduceRight",s,n)},shift(){return Qa(this,"shift")},some(s,n){return Vn(this,"some",s,n,void 0,arguments)},splice(...s){return Qa(this,"splice",s)},toReversed(){return Ma(this).toReversed()},toSorted(s){return Ma(this).toSorted(s)},toSpliced(...s){return Ma(this).toSpliced(...s)},unshift(...s){return Qa(this,"unshift",s)},values(){return Vt(this,"values",Qs)}};function Vt(s,n,a){const e=xt(s),t=e[n]();return e!==s&&!Pn(s)&&(t._next=t.next,t.next=()=>{const l=t._next();return l.done||(l.value=a(l.value)),l}),t}const Fu=Array.prototype;function Vn(s,n,a,e,t,l){const o=xt(s),p=o!==s&&!Pn(s),c=o[n];if(c!==Fu[n]){const i=c.apply(s,l);return p?Qs(i):i}let u=a;o!==s&&(p?u=function(i,d){return a.call(this,Qs(i),d,s)}:a.length>2&&(u=function(i,d){return a.call(this,i,d,s)}));const r=c.call(o,u,e);return p&&t?t(r):r}function Lo(s,n,a,e){const t=xt(s);let l=a;return t!==s&&(Pn(s)?a.length>3&&(l=function(o,p,c){return a.call(this,o,p,c,s)}):l=function(o,p,c){return a.call(this,o,Qs(p),c,s)}),t[n](l,...e)}function Ut(s,n,a){const e=xs(s);ln(e,"iterate",ge);const t=e[n](...a);return(t===-1||t===!1)&&Jl(a[0])?(a[0]=xs(a[0]),e[n](...a)):t}function Qa(s,n,a=[]){Zn(),Gl();const e=xs(s)[n].apply(s,a);return Wl(),sa(),e}const $u=Bl("__proto__,__v_isRef,__isVue"),Hc=new Set(Object.getOwnPropertyNames(Symbol).filter(s=>s!=="arguments"&&s!=="caller").map(s=>Symbol[s]).filter(fa));function Bu(s){fa(s)||(s=String(s));const n=xs(this);return ln(n,"has",s),n.hasOwnProperty(s)}class Gc{constructor(n=!1,a=!1){this._isReadonly=n,this._isShallow=a}get(n,a,e){if(a==="__v_skip")return n.__v_skip;const t=this._isReadonly,l=this._isShallow;if(a==="__v_isReactive")return!t;if(a==="__v_isReadonly")return t;if(a==="__v_isShallow")return l;if(a==="__v_raw")return e===(t?l?Yu:Yc:l?Xc:Kc).get(n)||Object.getPrototypeOf(n)===Object.getPrototypeOf(e)?n:void 0;const o=ds(n);if(!t){let c;if(o&&(c=Nu[a]))return c;if(a==="hasOwnProperty")return Bu}const p=Reflect.get(n,a,Fs(n)?n:e);if((fa(a)?Hc.has(a):$u(a))||(t||ln(n,"get",a),l))return p;if(Fs(p)){const c=o&&Vl(a)?p:p.value;return t&&Ds(c)?dl(c):c}return Ds(p)?t?dl(p):Se(p):p}}class Wc extends Gc{constructor(n=!1){super(!1,n)}set(n,a,e,t){let l=n[a];if(!this._isShallow){const c=ma(l);if(!Pn(e)&&!ma(e)&&(l=xs(l),e=xs(e)),!ds(n)&&Fs(l)&&!Fs(e))return c||(l.value=e),!0}const o=ds(n)&&Vl(a)?Number(a)<n.length:Es(n,a),p=Reflect.set(n,a,e,Fs(n)?n:t);return n===xs(t)&&(o?ha(e,l)&&Xn(n,"set",a,e):Xn(n,"add",a,e)),p}deleteProperty(n,a){const e=Es(n,a);n[a];const t=Reflect.deleteProperty(n,a);return t&&e&&Xn(n,"delete",a,void 0),t}has(n,a){const e=Reflect.has(n,a);return(!fa(a)||!Hc.has(a))&&ln(n,"has",a),e}ownKeys(n){return ln(n,"iterate",ds(n)?"length":Ta),Reflect.ownKeys(n)}}class qu extends Gc{constructor(n=!1){super(!0,n)}set(n,a){return!0}deleteProperty(n,a){return!0}}const zu=new Wc,Vu=new qu,Uu=new Wc(!0);const hl=s=>s,$e=s=>Reflect.getPrototypeOf(s);function Hu(s,n,a){return function(...e){const t=this.__v_raw,l=xs(t),o=qa(l),p=s==="entries"||s===Symbol.iterator&&o,c=s==="keys"&&o,u=t[s](...e),r=a?hl:n?ot:Qs;return!n&&ln(l,"iterate",c?ul:Ta),{next(){const{value:i,done:d}=u.next();return d?{value:i,done:d}:{value:p?[r(i[0]),r(i[1])]:r(i),done:d}},[Symbol.iterator](){return this}}}}function Be(s){return function(...n){return s==="delete"?!1:s==="clear"?void 0:this}}function Gu(s,n){const a={get(t){const l=this.__v_raw,o=xs(l),p=xs(t);s||(ha(t,p)&&ln(o,"get",t),ln(o,"get",p));const{has:c}=$e(o),u=n?hl:s?ot:Qs;if(c.call(o,t))return u(l.get(t));if(c.call(o,p))return u(l.get(p));l!==o&&l.get(t)},get size(){const t=this.__v_raw;return!s&&ln(xs(t),"iterate",Ta),t.size},has(t){const l=this.__v_raw,o=xs(l),p=xs(t);return s||(ha(t,p)&&ln(o,"has",t),ln(o,"has",p)),t===p?l.has(t):l.has(t)||l.has(p)},forEach(t,l){const o=this,p=o.__v_raw,c=xs(p),u=n?hl:s?ot:Qs;return!s&&ln(c,"iterate",Ta),p.forEach((r,i)=>t.call(l,u(r),u(i),o))}};return Ks(a,s?{add:Be("add"),set:Be("set"),delete:Be("delete"),clear:Be("clear")}:{add(t){!n&&!Pn(t)&&!ma(t)&&(t=xs(t));const l=xs(this);return $e(l).has.call(l,t)||(l.add(t),Xn(l,"add",t,t)),this},set(t,l){!n&&!Pn(l)&&!ma(l)&&(l=xs(l));const o=xs(this),{has:p,get:c}=$e(o);let u=p.call(o,t);u||(t=xs(t),u=p.call(o,t));const r=c.call(o,t);return o.set(t,l),u?ha(l,r)&&Xn(o,"set",t,l):Xn(o,"add",t,l),this},delete(t){const l=xs(this),{has:o,get:p}=$e(l);let c=o.call(l,t);c||(t=xs(t),c=o.call(l,t)),p&&p.call(l,t);const u=l.delete(t);return c&&Xn(l,"delete",t,void 0),u},clear(){const t=xs(this),l=t.size!==0,o=t.clear();return l&&Xn(t,"clear",void 0,void 0),o}}),["keys","values","entries",Symbol.iterator].forEach(t=>{a[t]=Hu(t,s,n)}),a}function Yl(s,n){const a=Gu(s,n);return(e,t,l)=>t==="__v_isReactive"?!s:t==="__v_isReadonly"?s:t==="__v_raw"?e:Reflect.get(Es(a,t)&&t in e?a:e,t,l)}const Wu={get:Yl(!1,!1)},Ku={get:Yl(!1,!0)},Xu={get:Yl(!0,!1)};const Kc=new WeakMap,Xc=new WeakMap,Yc=new WeakMap,Yu=new WeakMap;function Qu(s){switch(s){case"Object":case"Array":return 1;case"Map":case"Set":case"WeakMap":case"WeakSet":return 2;default:return 0}}function Ju(s){return s.__v_skip||!Object.isExtensible(s)?0:Qu(Cu(s))}function Se(s){return ma(s)?s:Ql(s,!1,zu,Wu,Kc)}function Qc(s){return Ql(s,!1,Uu,Ku,Xc)}function dl(s){return Ql(s,!0,Vu,Xu,Yc)}function Ql(s,n,a,e,t){if(!Ds(s)||s.__v_raw&&!(n&&s.__v_isReactive))return s;const l=Ju(s);if(l===0)return s;const o=t.get(s);if(o)return o;const p=new Proxy(s,l===2?e:a);return t.set(s,p),p}function da(s){return ma(s)?da(s.__v_raw):!!(s&&s.__v_isReactive)}function ma(s){return!!(s&&s.__v_isReadonly)}function Pn(s){return!!(s&&s.__v_isShallow)}function Jl(s){return s?!!s.__v_raw:!1}function xs(s){const n=s&&s.__v_raw;return n?xs(n):s}function Zl(s){return!Es(s,"__v_skip")&&Object.isExtensible(s)&&Rc(s,"__v_skip",!0),s}const Qs=s=>Ds(s)?Se(s):s,ot=s=>Ds(s)?dl(s):s;function Fs(s){return s?s.__v_isRef===!0:!1}function os(s){return Jc(s,!1)}function so(s){return Jc(s,!0)}function Jc(s,n){return Fs(s)?s:new Zu(s,n)}class Zu{constructor(n,a){this.dep=new Xl,this.__v_isRef=!0,this.__v_isShallow=!1,this._rawValue=a?n:xs(n),this._value=a?n:Qs(n),this.__v_isShallow=a}get value(){return this.dep.track(),this._value}set value(n){const a=this._rawValue,e=this.__v_isShallow||Pn(n)||ma(n);n=e?n:xs(n),ha(n,a)&&(this._rawValue=n,this._value=e?n:Qs(n),this.dep.trigger())}}function cs(s){return Fs(s)?s.value:s}function sh(s){return js(s)?s():cs(s)}const nh={get:(s,n,a)=>n==="__v_raw"?s:cs(Reflect.get(s,n,a)),set:(s,n,a,e)=>{const t=s[n];return Fs(t)&&!Fs(a)?(t.value=a,!0):Reflect.set(s,n,a,e)}};function Zc(s){return da(s)?s:new Proxy(s,nh)}function ah(s){const n=ds(s)?new Array(s.length):{};for(const a in s)n[a]=th(s,a);return n}class eh{constructor(n,a,e){this._object=n,this._key=a,this._defaultValue=e,this.__v_isRef=!0,this._value=void 0}get value(){const n=this._object[this._key];return this._value=n===void 0?this._defaultValue:n}set value(n){this._object[this._key]=n}get dep(){return Iu(xs(this._object),this._key)}}function th(s,n,a){const e=s[n];return Fs(e)?e:new eh(s,n,a)}class lh{constructor(n,a,e){this.fn=n,this.setter=a,this._value=void 0,this.dep=new Xl(this),this.__v_isRef=!0,this.deps=void 0,this.depsTail=void 0,this.flags=16,this.globalVersion=fe-1,this.next=void 0,this.effect=this,this.__v_isReadonly=!a,this.isSSR=e}notify(){if(this.flags|=16,!(this.flags&8)&&Ls!==this)return $c(this,!0),!0}get value(){const n=this.dep.track();return zc(this),n&&(n.version=this.dep.version),this._value}set value(n){this.setter&&this.setter(n)}}function oh(s,n,a=!1){let e,t;return js(s)?e=s:(e=s.get,t=s.set),new lh(e,t,a)}const qe={},pt=new WeakMap;let Pa;function ph(s,n=!1,a=Pa){if(a){let e=pt.get(a);e||pt.set(a,e=[]),e.push(s)}}function ch(s,n,a=Ss){const{immediate:e,deep:t,once:l,scheduler:o,augmentJob:p,call:c}=a,u=f=>t?f:Pn(f)||t===!1||t===0?Yn(f,1):Yn(f);let r,i,d,h,y=!1,j=!1;if(Fs(s)?(i=()=>s.value,y=Pn(s)):da(s)?(i=()=>u(s),y=!0):ds(s)?(j=!0,y=s.some(f=>da(f)||Pn(f)),i=()=>s.map(f=>{if(Fs(f))return f.value;if(da(f))return u(f);if(js(f))return c?c(f,2):f()})):js(s)?n?i=c?()=>c(s,2):s:i=()=>{if(d){Zn();try{d()}finally{sa()}}const f=Pa;Pa=r;try{return c?c(s,3,[h]):s(h)}finally{Pa=f}}:i=$n,n&&t){const f=i,x=t===!0?1/0:t;i=()=>Yn(f(),x)}const T=Ic(),g=()=>{r.stop(),T&&T.active&&zl(T.effects,r)};if(l&&n){const f=n;n=(...x)=>{f(...x),g()}}let m=j?new Array(s.length).fill(qe):qe;const w=f=>{if(!(!(r.flags&1)||!r.dirty&&!f))if(n){const x=r.run();if(t||y||(j?x.some((k,_)=>ha(k,m[_])):ha(x,m))){d&&d();const k=Pa;Pa=r;try{const _=[x,m===qe?void 0:j&&m[0]===qe?[]:m,h];m=x,c?c(n,3,_):n(..._)}finally{Pa=k}}}else r.run()};return p&&p(w),r=new Nc(i),r.scheduler=o?()=>o(w,!1):w,h=f=>ph(f,!1,r),d=r.onStop=()=>{const f=pt.get(r);if(f){if(c)c(f,4);else for(const x of f)x();pt.delete(r)}},n?e?w(!0):m=r.run():o?o(w.bind(null,!0),!0):r.run(),g.pause=r.pause.bind(r),g.resume=r.resume.bind(r),g.stop=g,g}function Yn(s,n=1/0,a){if(n<=0||!Ds(s)||s.__v_skip||(a=a||new Map,(a.get(s)||0)>=n))return s;if(a.set(s,n),n--,Fs(s))Yn(s.value,n,a);else if(ds(s))for(let e=0;e<s.length;e++)Yn(s[e],n,a);else if(Ec(s)||qa(s))s.forEach(e=>{Yn(e,n,a)});else if(Ac(s)){for(const e in s)Yn(s[e],n,a);for(const e of Object.getOwnPropertySymbols(s))Object.prototype.propertyIsEnumerable.call(s,e)&&Yn(s[e],n,a)}return s}/**
* @vue/runtime-core v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/function Te(s,n,a,e){try{return e?s(...e):s()}catch(t){Ae(t,n,a)}}function Ln(s,n,a,e){if(js(s)){const t=Te(s,n,a,e);return t&&Sc(t)&&t.catch(l=>{Ae(l,n,a)}),t}if(ds(s)){const t=[];for(let l=0;l<s.length;l++)t.push(Ln(s[l],n,a,e));return t}}function Ae(s,n,a,e=!0){const t=n?n.vnode:null,{errorHandler:l,throwUnhandledErrorInProduction:o}=n&&n.appContext.config||Ss;if(n){let p=n.parent;const c=n.proxy,u=`https://vuejs.org/error-reference/#runtime-${a}`;for(;p;){const r=p.ec;if(r){for(let i=0;i<r.length;i++)if(r[i](s,c,u)===!1)return}p=p.parent}if(l){Zn(),Te(l,null,10,[s,c,u]),sa();return}}rh(s,a,t,e,o)}function rh(s,n,a,e=!0,t=!1){if(t)throw s;console.error(s)}const hn=[];let Nn=-1;const za=[];let ca=null,Ia=0;const sr=Promise.resolve();let ct=null;function Js(s){const n=ct||sr;return s?n.then(this?s.bind(this):s):n}function ih(s){let n=Nn+1,a=hn.length;for(;n<a;){const e=n+a>>>1,t=hn[e],l=be(t);l<s||l===s&&t.flags&2?n=e+1:a=e}return n}function no(s){if(!(s.flags&1)){const n=be(s),a=hn[hn.length-1];!a||!(s.flags&2)&&n>=be(a)?hn.push(s):hn.splice(ih(n),0,s),s.flags|=1,nr()}}function nr(){ct||(ct=sr.then(er))}function uh(s){ds(s)?za.push(...s):ca&&s.id===-1?ca.splice(Ia+1,0,s):s.flags&1||(za.push(s),s.flags|=1),nr()}function Mo(s,n,a=Nn+1){for(;a<hn.length;a++){const e=hn[a];if(e&&e.flags&2){if(s&&e.id!==s.uid)continue;hn.splice(a,1),a--,e.flags&4&&(e.flags&=-2),e(),e.flags&4||(e.flags&=-2)}}}function ar(s){if(za.length){const n=[...new Set(za)].sort((a,e)=>be(a)-be(e));if(za.length=0,ca){ca.push(...n);return}for(ca=n,Ia=0;Ia<ca.length;Ia++){const a=ca[Ia];a.flags&4&&(a.flags&=-2),a.flags&8||a(),a.flags&=-2}ca=null,Ia=0}}const be=s=>s.id==null?s.flags&2?-1:1/0:s.id;function er(s){try{for(Nn=0;Nn<hn.length;Nn++){const n=hn[Nn];n&&!(n.flags&8)&&(n.flags&4&&(n.flags&=-2),Te(n,n.i,n.i?15:14),n.flags&4||(n.flags&=-2))}}finally{for(;Nn<hn.length;Nn++){const n=hn[Nn];n&&(n.flags&=-2)}Nn=-1,hn.length=0,ar(),ct=null,(hn.length||za.length)&&er()}}let gn=null,tr=null;function rt(s){const n=gn;return gn=s,tr=s&&s.type.__scopeId||null,n}function vn(s,n=gn,a){if(!n||s._n)return s;const e=(...t)=>{e._d&&ht(-1);const l=rt(n);let o;try{o=s(...t)}finally{rt(l),e._d&&ht(1)}return o};return e._n=!0,e._c=!0,e._d=!0,e}function Ua(s,n){if(gn===null)return s;const a=St(gn),e=s.dirs||(s.dirs=[]);for(let t=0;t<n.length;t++){let[l,o,p,c=Ss]=n[t];l&&(js(l)&&(l={mounted:l,updated:l}),l.deep&&Yn(o),e.push({dir:l,instance:a,value:o,oldValue:void 0,arg:p,modifiers:c}))}return s}function wa(s,n,a,e){const t=s.dirs,l=n&&n.dirs;for(let o=0;o<t.length;o++){const p=t[o];l&&(p.oldValue=l[o].value);let c=p.dir[e];c&&(Zn(),Ln(c,a,8,[s.el,p,s,n]),sa())}}const lr=Symbol("_vte"),or=s=>s.__isTeleport,pe=s=>s&&(s.disabled||s.disabled===""),Do=s=>s&&(s.defer||s.defer===""),Oo=s=>typeof SVGElement<"u"&&s instanceof SVGElement,Io=s=>typeof MathMLElement=="function"&&s instanceof MathMLElement,ml=(s,n)=>{const a=s&&s.to;return $s(a)?n?n(a):null:a},pr={name:"Teleport",__isTeleport:!0,process(s,n,a,e,t,l,o,p,c,u){const{mc:r,pc:i,pbc:d,o:{insert:h,querySelector:y,createText:j,createComment:T}}=u,g=pe(n.props);let{shapeFlag:m,children:w,dynamicChildren:f}=n;if(s==null){const x=n.el=j(""),k=n.anchor=j("");h(x,a,e),h(k,a,e);const _=(S,I)=>{m&16&&r(w,S,I,t,l,o,p,c)},D=()=>{const S=n.target=ml(n.props,y),I=cr(S,n,j,h);S&&(o!=="svg"&&Oo(S)?o="svg":o!=="mathml"&&Io(S)&&(o="mathml"),t&&t.isCE&&(t.ce._teleportTargets||(t.ce._teleportTargets=new Set)).add(S),g||(_(S,I),Qe(n,!1)))};g&&(_(a,k),Qe(n,!0)),Do(n.props)?(n.el.__isMounted=!1,rn(()=>{D(),delete n.el.__isMounted},l)):D()}else{if(Do(n.props)&&s.el.__isMounted===!1){rn(()=>{pr.process(s,n,a,e,t,l,o,p,c,u)},l);return}n.el=s.el,n.targetStart=s.targetStart;const x=n.anchor=s.anchor,k=n.target=s.target,_=n.targetAnchor=s.targetAnchor,D=pe(s.props),S=D?a:k,I=D?x:_;if(o==="svg"||Oo(k)?o="svg":(o==="mathml"||Io(k))&&(o="mathml"),f?(d(s.dynamicChildren,f,S,t,l,o,p),ro(s,n,!0)):c||i(s,n,S,I,t,l,o,p,!1),g)D?n.props&&s.props&&n.props.to!==s.props.to&&(n.props.to=s.props.to):ze(n,a,x,u,1);else if((n.props&&n.props.to)!==(s.props&&s.props.to)){const J=n.target=ml(n.props,y);J&&ze(n,J,null,u,0)}else D&&ze(n,k,_,u,1);Qe(n,g)}},remove(s,n,a,{um:e,o:{remove:t}},l){const{shapeFlag:o,children:p,anchor:c,targetStart:u,targetAnchor:r,target:i,props:d}=s;if(i&&(t(u),t(r)),l&&t(c),o&16){const h=l||!pe(d);for(let y=0;y<p.length;y++){const j=p[y];e(j,n,a,h,!!j.dynamicChildren)}}},move:ze,hydrate:hh};function ze(s,n,a,{o:{insert:e},m:t},l=2){l===0&&e(s.targetAnchor,n,a);const{el:o,anchor:p,shapeFlag:c,children:u,props:r}=s,i=l===2;if(i&&e(o,n,a),(!i||pe(r))&&c&16)for(let d=0;d<u.length;d++)t(u[d],n,a,2);i&&e(p,n,a)}function hh(s,n,a,e,t,l,{o:{nextSibling:o,parentNode:p,querySelector:c,insert:u,createText:r}},i){function d(j,T,g,m){T.anchor=i(o(j),T,p(j),a,e,t,l),T.targetStart=g,T.targetAnchor=m}const h=n.target=ml(n.props,c),y=pe(n.props);if(h){const j=h._lpa||h.firstChild;if(n.shapeFlag&16)if(y)d(s,n,j,j&&o(j));else{n.anchor=o(s);let T=j;for(;T;){if(T&&T.nodeType===8){if(T.data==="teleport start anchor")n.targetStart=T;else if(T.data==="teleport anchor"){n.targetAnchor=T,h._lpa=n.targetAnchor&&o(n.targetAnchor);break}}T=o(T)}n.targetAnchor||cr(h,n,r,u),i(j&&o(j),n,h,a,e,t,l)}Qe(n,y)}else y&&n.shapeFlag&16&&d(s,n,s,o(s));return n.anchor&&o(n.anchor)}const Cw=pr;function Qe(s,n){const a=s.ctx;if(a&&a.ut){let e,t;for(n?(e=s.el,t=s.anchor):(e=s.targetStart,t=s.targetAnchor);e&&e!==t;)e.nodeType===1&&e.setAttribute("data-v-owner",a.uid),e=e.nextSibling;a.ut()}}function cr(s,n,a,e){const t=n.targetStart=a(""),l=n.targetAnchor=a("");return t[lr]=l,s&&(e(t,s),e(l,s)),l}const Kn=Symbol("_leaveCb"),Ve=Symbol("_enterCb");function dh(){const s={isMounted:!1,isLeaving:!1,isUnmounting:!1,leavingVNodes:new Map};return wn(()=>{s.isMounted=!0}),qn(()=>{s.isUnmounting=!0}),s}const kn=[Function,Array],rr={mode:String,appear:Boolean,persisted:Boolean,onBeforeEnter:kn,onEnter:kn,onAfterEnter:kn,onEnterCancelled:kn,onBeforeLeave:kn,onLeave:kn,onAfterLeave:kn,onLeaveCancelled:kn,onBeforeAppear:kn,onAppear:kn,onAfterAppear:kn,onAppearCancelled:kn},ir=s=>{const n=s.subTree;return n.component?ir(n.component):n},mh={name:"BaseTransition",props:rr,setup(s,{slots:n}){const a=ea(),e=dh();return()=>{const t=n.default&&dr(n.default(),!0);if(!t||!t.length)return;const l=ur(t),o=xs(s),{mode:p}=o;if(e.isLeaving)return Ht(l);const c=No(l);if(!c)return Ht(l);let u=jl(c,o,e,a,i=>u=i);c.type!==dn&&_e(c,u);let r=a.subTree&&No(a.subTree);if(r&&r.type!==dn&&!Ea(r,c)&&ir(a).type!==dn){let i=jl(r,o,e,a);if(_e(r,i),p==="out-in"&&c.type!==dn)return e.isLeaving=!0,i.afterLeave=()=>{e.isLeaving=!1,a.job.flags&8||a.update(),delete i.afterLeave,r=void 0},Ht(l);p==="in-out"&&c.type!==dn?i.delayLeave=(d,h,y)=>{const j=hr(e,r);j[String(r.key)]=r,d[Kn]=()=>{h(),d[Kn]=void 0,delete u.delayedLeave,r=void 0},u.delayedLeave=()=>{y(),delete u.delayedLeave,r=void 0}}:r=void 0}else r&&(r=void 0);return l}}};function ur(s){let n=s[0];if(s.length>1){for(const a of s)if(a.type!==dn){n=a;break}}return n}const jh=mh;function hr(s,n){const{leavingVNodes:a}=s;let e=a.get(n.type);return e||(e=Object.create(null),a.set(n.type,e)),e}function jl(s,n,a,e,t){const{appear:l,mode:o,persisted:p=!1,onBeforeEnter:c,onEnter:u,onAfterEnter:r,onEnterCancelled:i,onBeforeLeave:d,onLeave:h,onAfterLeave:y,onLeaveCancelled:j,onBeforeAppear:T,onAppear:g,onAfterAppear:m,onAppearCancelled:w}=n,f=String(s.key),x=hr(a,s),k=(S,I)=>{S&&Ln(S,e,9,I)},_=(S,I)=>{const J=I[1];k(S,I),ds(S)?S.every(U=>U.length<=1)&&J():S.length<=1&&J()},D={mode:o,persisted:p,beforeEnter(S){let I=c;if(!a.isMounted)if(l)I=T||c;else return;S[Kn]&&S[Kn](!0);const J=x[f];J&&Ea(s,J)&&J.el[Kn]&&J.el[Kn](),k(I,[S])},enter(S){let I=u,J=r,U=i;if(!a.isMounted)if(l)I=g||u,J=m||r,U=w||i;else return;let F=!1;const G=S[Ve]=is=>{F||(F=!0,is?k(U,[S]):k(J,[S]),D.delayedLeave&&D.delayedLeave(),S[Ve]=void 0)};I?_(I,[S,G]):G()},leave(S,I){const J=String(s.key);if(S[Ve]&&S[Ve](!0),a.isUnmounting)return I();k(d,[S]);let U=!1;const F=S[Kn]=G=>{U||(U=!0,I(),G?k(j,[S]):k(y,[S]),S[Kn]=void 0,x[J]===s&&delete x[J])};x[J]=s,h?_(h,[S,F]):F()},clone(S){const I=jl(S,n,a,e,t);return t&&t(I),I}};return D}function Ht(s){if(Re(s))return s=ja(s),s.children=null,s}function No(s){if(!Re(s))return or(s.type)&&s.children?ur(s.children):s;if(s.component)return s.component.subTree;const{shapeFlag:n,children:a}=s;if(a){if(n&16)return a[0];if(n&32&&js(a.default))return a.default()}}function _e(s,n){s.shapeFlag&6&&s.component?(s.transition=n,_e(s.component.subTree,n)):s.shapeFlag&128?(s.ssContent.transition=n.clone(s.ssContent),s.ssFallback.transition=n.clone(s.ssFallback)):s.transition=n}function dr(s,n=!1,a){let e=[],t=0;for(let l=0;l<s.length;l++){let o=s[l];const p=a==null?o.key:String(a)+String(o.key!=null?o.key:l);o.type===bs?(o.patchFlag&128&&t++,e=e.concat(dr(o.children,n,p))):(n||o.type!==dn)&&e.push(p!=null?ja(o,{key:p}):o)}if(t>1)for(let l=0;l<e.length;l++)e[l].patchFlag=-2;return e}function Gs(s,n){return js(s)?Ks({name:s.name},n,{setup:s}):s}function ao(s){s.ids=[s.ids[0]+s.ids[2]+++"-",0,0]}function Je(s){const n=ea(),a=so(null);if(n){const t=n.refs===Ss?n.refs={}:n.refs;Object.defineProperty(t,s,{enumerable:!0,get:()=>a.value,set:l=>a.value=l})}return a}const it=new WeakMap;function ce(s,n,a,e,t=!1){if(ds(s)){s.forEach((y,j)=>ce(y,n&&(ds(n)?n[j]:n),a,e,t));return}if(re(e)&&!t){e.shapeFlag&512&&e.type.__asyncResolved&&e.component.subTree.component&&ce(s,n,a,e.component.subTree);return}const l=e.shapeFlag&4?St(e.component):e.el,o=t?null:l,{i:p,r:c}=s,u=n&&n.r,r=p.refs===Ss?p.refs={}:p.refs,i=p.setupState,d=xs(i),h=i===Ss?Pc:y=>Es(d,y);if(u!=null&&u!==c){if(Fo(n),$s(u))r[u]=null,h(u)&&(i[u]=null);else if(Fs(u)){u.value=null;const y=n;y.k&&(r[y.k]=null)}}if(js(c))Te(c,p,12,[o,r]);else{const y=$s(c),j=Fs(c);if(y||j){const T=()=>{if(s.f){const g=y?h(c)?i[c]:r[c]:c.value;if(t)ds(g)&&zl(g,l);else if(ds(g))g.includes(l)||g.push(l);else if(y)r[c]=[l],h(c)&&(i[c]=r[c]);else{const m=[l];c.value=m,s.k&&(r[s.k]=m)}}else y?(r[c]=o,h(c)&&(i[c]=o)):j&&(c.value=o,s.k&&(r[s.k]=o))};if(o){const g=()=>{T(),it.delete(s)};g.id=-1,it.set(s,g),rn(g,a)}else Fo(s),T()}}}function Fo(s){const n=it.get(s);n&&(n.flags|=8,it.delete(s))}const $o=s=>s.nodeType===8;kt().requestIdleCallback;kt().cancelIdleCallback;function fh(s,n){if($o(s)&&s.data==="["){let a=1,e=s.nextSibling;for(;e;){if(e.nodeType===1){if(n(e)===!1)break}else if($o(e))if(e.data==="]"){if(--a===0)break}else e.data==="["&&a++;e=e.nextSibling}}else n(s)}const re=s=>!!s.type.__asyncLoader;function gh(s){js(s)&&(s={loader:s});const{loader:n,loadingComponent:a,errorComponent:e,delay:t=200,hydrate:l,timeout:o,suspensible:p=!0,onError:c}=s;let u=null,r,i=0;const d=()=>(i++,u=null,h()),h=()=>{let y;return u||(y=u=n().catch(j=>{if(j=j instanceof Error?j:new Error(String(j)),c)return new Promise((T,g)=>{c(j,()=>T(d()),()=>g(j),i+1)});throw j}).then(j=>y!==u&&u?u:(j&&(j.__esModule||j[Symbol.toStringTag]==="Module")&&(j=j.default),r=j,j)))};return Gs({name:"AsyncComponentWrapper",__asyncLoader:h,__asyncHydrate(y,j,T){let g=!1;(j.bu||(j.bu=[])).push(()=>g=!0);const m=()=>{g||T()},w=l?()=>{const f=l(m,x=>fh(y,x));f&&(j.bum||(j.bum=[])).push(f)}:m;r?w():h().then(()=>!j.isUnmounted&&w())},get __asyncResolved(){return r},setup(){const y=Zs;if(ao(y),r)return()=>Ue(r,y);const j=w=>{u=null,Ae(w,y,13,!e)};if(p&&y.suspense||Ha)return h().then(w=>()=>Ue(w,y)).catch(w=>(j(w),()=>e?ws(e,{error:w}):null));const T=os(!1),g=os(),m=os(!!t);return t&&setTimeout(()=>{m.value=!1},t),o!=null&&setTimeout(()=>{if(!T.value&&!g.value){const w=new Error(`Async component timed out after ${o}ms.`);j(w),g.value=w}},o),h().then(()=>{T.value=!0,y.parent&&Re(y.parent.vnode)&&y.parent.update()}).catch(w=>{j(w),g.value=w}),()=>{if(T.value&&r)return Ue(r,y);if(g.value&&e)return ws(e,{error:g.value});if(a&&!m.value)return Ue(a,y)}}})}function Ue(s,n){const{ref:a,props:e,children:t,ce:l}=n.vnode,o=ws(s,e,t);return o.ref=a,o.ce=l,delete n.vnode.ce,o}const Re=s=>s.type.__isKeepAlive;function mr(s,n){fr(s,"a",n)}function jr(s,n){fr(s,"da",n)}function fr(s,n,a=Zs){const e=s.__wdc||(s.__wdc=()=>{let t=a;for(;t;){if(t.isDeactivated)return;t=t.parent}return s()});if(Pt(n,e,a),a){let t=a.parent;for(;t&&t.parent;)Re(t.parent.vnode)&&bh(e,n,a,t),t=t.parent}}function bh(s,n,a,e){const t=Pt(n,s,e,!0);eo(()=>{zl(e[n],t)},a)}function Pt(s,n,a=Zs,e=!1){if(a){const t=a[s]||(a[s]=[]),l=n.__weh||(n.__weh=(...o)=>{Zn();const p=De(a),c=Ln(n,a,s,o);return p(),sa(),c});return e?t.unshift(l):t.push(l),l}}const aa=s=>(n,a=Zs)=>{(!Ha||s==="sp")&&Pt(s,(...e)=>n(...e),a)},_h=aa("bm"),wn=aa("m"),yh=aa("bu"),vh=aa("u"),qn=aa("bum"),eo=aa("um"),wh=aa("sp"),Ch=aa("rtg"),kh=aa("rtc");function xh(s,n=Zs){Pt("ec",s,n)}const to="components",Ph="directives";function Le(s,n){return oo(to,s,!0,n)||s}const gr=Symbol.for("v-ndc");function Eh(s){return $s(s)?oo(to,s,!1)||s:s||gr}function lo(s){return oo(Ph,s)}function oo(s,n,a=!0,e=!1){const t=gn||Zs;if(t){const l=t.type;if(s===to){const p=jd(l,!1);if(p&&(p===n||p===Tn(n)||p===Ct(Tn(n))))return l}const o=Bo(t[s]||l[s],n)||Bo(t.appContext[s],n);return!o&&e?l:o}}function Bo(s,n){return s&&(s[n]||s[Tn(n)]||s[Ct(Tn(n))])}function Ms(s,n,a,e){let t;const l=a,o=ds(s);if(o||$s(s)){const p=o&&da(s);let c=!1,u=!1;p&&(c=!Pn(s),u=ma(s),s=xt(s)),t=new Array(s.length);for(let r=0,i=s.length;r<i;r++)t[r]=n(c?u?ot(Qs(s[r])):Qs(s[r]):s[r],r,void 0,l)}else if(typeof s=="number"){t=new Array(s);for(let p=0;p<s;p++)t[p]=n(p+1,p,void 0,l)}else if(Ds(s))if(s[Symbol.iterator])t=Array.from(s,(p,c)=>n(p,c,void 0,l));else{const p=Object.keys(s);t=new Array(p.length);for(let c=0,u=p.length;c<u;c++){const r=p[c];t[c]=n(s[r],r,c,l)}}else t=[];return t}const fl=s=>s?Ir(s)?St(s):fl(s.parent):null,ie=Ks(Object.create(null),{$:s=>s,$el:s=>s.vnode.el,$data:s=>s.data,$props:s=>s.props,$attrs:s=>s.attrs,$slots:s=>s.slots,$refs:s=>s.refs,$parent:s=>fl(s.parent),$root:s=>fl(s.root),$host:s=>s.ce,$emit:s=>s.emit,$options:s=>_r(s),$forceUpdate:s=>s.f||(s.f=()=>{no(s.update)}),$nextTick:s=>s.n||(s.n=Js.bind(s.proxy)),$watch:s=>Xh.bind(s)}),Gt=(s,n)=>s!==Ss&&!s.__isScriptSetup&&Es(s,n),Sh={get({_:s},n){if(n==="__v_skip")return!0;const{ctx:a,setupState:e,data:t,props:l,accessCache:o,type:p,appContext:c}=s;let u;if(n[0]!=="$"){const h=o[n];if(h!==void 0)switch(h){case 1:return e[n];case 2:return t[n];case 4:return a[n];case 3:return l[n]}else{if(Gt(e,n))return o[n]=1,e[n];if(t!==Ss&&Es(t,n))return o[n]=2,t[n];if((u=s.propsOptions[0])&&Es(u,n))return o[n]=3,l[n];if(a!==Ss&&Es(a,n))return o[n]=4,a[n];gl&&(o[n]=0)}}const r=ie[n];let i,d;if(r)return n==="$attrs"&&ln(s.attrs,"get",""),r(s);if((i=p.__cssModules)&&(i=i[n]))return i;if(a!==Ss&&Es(a,n))return o[n]=4,a[n];if(d=c.config.globalProperties,Es(d,n))return d[n]},set({_:s},n,a){const{data:e,setupState:t,ctx:l}=s;return Gt(t,n)?(t[n]=a,!0):e!==Ss&&Es(e,n)?(e[n]=a,!0):Es(s.props,n)||n[0]==="$"&&n.slice(1)in s?!1:(l[n]=a,!0)},has({_:{data:s,setupState:n,accessCache:a,ctx:e,appContext:t,propsOptions:l,type:o}},p){let c,u;return!!(a[p]||s!==Ss&&p[0]!=="$"&&Es(s,p)||Gt(n,p)||(c=l[0])&&Es(c,p)||Es(e,p)||Es(ie,p)||Es(t.config.globalProperties,p)||(u=o.__cssModules)&&u[p])},defineProperty(s,n,a){return a.get!=null?s._.accessCache[n]=0:Es(a,"value")&&this.set(s,n,a.value,null),Reflect.defineProperty(s,n,a)}};function qo(s){return ds(s)?s.reduce((n,a)=>(n[a]=null,n),{}):s}let gl=!0;function Th(s){const n=_r(s),a=s.proxy,e=s.ctx;gl=!1,n.beforeCreate&&zo(n.beforeCreate,s,"bc");const{data:t,computed:l,methods:o,watch:p,provide:c,inject:u,created:r,beforeMount:i,mounted:d,beforeUpdate:h,updated:y,activated:j,deactivated:T,beforeDestroy:g,beforeUnmount:m,destroyed:w,unmounted:f,render:x,renderTracked:k,renderTriggered:_,errorCaptured:D,serverPrefetch:S,expose:I,inheritAttrs:J,components:U,directives:F,filters:G}=n;if(u&&Ah(u,e,null),o)for(const Y in o){const us=o[Y];js(us)&&(e[Y]=us.bind(a))}if(t){const Y=t.call(a,a);Ds(Y)&&(s.data=Se(Y))}if(gl=!0,l)for(const Y in l){const us=l[Y],Ns=js(us)?us.bind(a,a):js(us.get)?us.get.bind(a,a):$n,Xs=!js(us)&&js(us.set)?us.set.bind(a):$n,Rs=ss({get:Ns,set:Xs});Object.defineProperty(e,Y,{enumerable:!0,configurable:!0,get:()=>Rs.value,set:zs=>Rs.value=zs})}if(p)for(const Y in p)br(p[Y],e,a,Y);if(c){const Y=js(c)?c.call(a):c;Reflect.ownKeys(Y).forEach(us=>{Ze(us,Y[us])})}r&&zo(r,s,"c");function ps(Y,us){ds(us)?us.forEach(Ns=>Y(Ns.bind(a))):us&&Y(us.bind(a))}if(ps(_h,i),ps(wn,d),ps(yh,h),ps(vh,y),ps(mr,j),ps(jr,T),ps(xh,D),ps(kh,k),ps(Ch,_),ps(qn,m),ps(eo,f),ps(wh,S),ds(I))if(I.length){const Y=s.exposed||(s.exposed={});I.forEach(us=>{Object.defineProperty(Y,us,{get:()=>a[us],set:Ns=>a[us]=Ns,enumerable:!0})})}else s.exposed||(s.exposed={});x&&s.render===$n&&(s.render=x),J!=null&&(s.inheritAttrs=J),U&&(s.components=U),F&&(s.directives=F),S&&ao(s)}function Ah(s,n,a=$n){ds(s)&&(s=bl(s));for(const e in s){const t=s[e];let l;Ds(t)?"default"in t?l=jn(t.from||e,t.default,!0):l=jn(t.from||e):l=jn(t),Fs(l)?Object.defineProperty(n,e,{enumerable:!0,configurable:!0,get:()=>l.value,set:o=>l.value=o}):n[e]=l}}function zo(s,n,a){Ln(ds(s)?s.map(e=>e.bind(n.proxy)):s.bind(n.proxy),n,a)}function br(s,n,a,e){let t=e.includes(".")?Rr(a,e):()=>a[e];if($s(s)){const l=n[s];js(l)&&Hs(t,l)}else if(js(s))Hs(t,s.bind(a));else if(Ds(s))if(ds(s))s.forEach(l=>br(l,n,a,e));else{const l=js(s.handler)?s.handler.bind(a):n[s.handler];js(l)&&Hs(t,l,s)}}function _r(s){const n=s.type,{mixins:a,extends:e}=n,{mixins:t,optionsCache:l,config:{optionMergeStrategies:o}}=s.appContext,p=l.get(n);let c;return p?c=p:!t.length&&!a&&!e?c=n:(c={},t.length&&t.forEach(u=>ut(c,u,o,!0)),ut(c,n,o)),Ds(n)&&l.set(n,c),c}function ut(s,n,a,e=!1){const{mixins:t,extends:l}=n;l&&ut(s,l,a,!0),t&&t.forEach(o=>ut(s,o,a,!0));for(const o in n)if(!(e&&o==="expose")){const p=Rh[o]||a&&a[o];s[o]=p?p(s[o],n[o]):n[o]}return s}const Rh={data:Vo,props:Uo,emits:Uo,methods:ee,computed:ee,beforeCreate:pn,created:pn,beforeMount:pn,mounted:pn,beforeUpdate:pn,updated:pn,beforeDestroy:pn,beforeUnmount:pn,destroyed:pn,unmounted:pn,activated:pn,deactivated:pn,errorCaptured:pn,serverPrefetch:pn,components:ee,directives:ee,watch:Mh,provide:Vo,inject:Lh};function Vo(s,n){return n?s?function(){return Ks(js(s)?s.call(this,this):s,js(n)?n.call(this,this):n)}:n:s}function Lh(s,n){return ee(bl(s),bl(n))}function bl(s){if(ds(s)){const n={};for(let a=0;a<s.length;a++)n[s[a]]=s[a];return n}return s}function pn(s,n){return s?[...new Set([].concat(s,n))]:n}function ee(s,n){return s?Ks(Object.create(null),s,n):n}function Uo(s,n){return s?ds(s)&&ds(n)?[...new Set([...s,...n])]:Ks(Object.create(null),qo(s),qo(n??{})):n}function Mh(s,n){if(!s)return n;if(!n)return s;const a=Ks(Object.create(null),s);for(const e in n)a[e]=pn(s[e],n[e]);return a}function yr(){return{app:null,config:{isNativeTag:Pc,performance:!1,globalProperties:{},optionMergeStrategies:{},errorHandler:void 0,warnHandler:void 0,compilerOptions:{}},mixins:[],components:{},directives:{},provides:Object.create(null),optionsCache:new WeakMap,propsCache:new WeakMap,emitsCache:new WeakMap}}let Dh=0;function Oh(s,n){return function(e,t=null){js(e)||(e=Ks({},e)),t!=null&&!Ds(t)&&(t=null);const l=yr(),o=new WeakSet,p=[];let c=!1;const u=l.app={_uid:Dh++,_component:e,_props:t,_container:null,_context:l,_instance:null,version:gd,get config(){return l.config},set config(r){},use(r,...i){return o.has(r)||(r&&js(r.install)?(o.add(r),r.install(u,...i)):js(r)&&(o.add(r),r(u,...i))),u},mixin(r){return l.mixins.includes(r)||l.mixins.push(r),u},component(r,i){return i?(l.components[r]=i,u):l.components[r]},directive(r,i){return i?(l.directives[r]=i,u):l.directives[r]},mount(r,i,d){if(!c){const h=u._ceVNode||ws(e,t);return h.appContext=l,d===!0?d="svg":d===!1&&(d=void 0),s(h,r,d),c=!0,u._container=r,r.__vue_app__=u,St(h.component)}},onUnmount(r){p.push(r)},unmount(){c&&(Ln(p,u._instance,16),s(null,u._container),delete u._container.__vue_app__)},provide(r,i){return l.provides[r]=i,u},runWithContext(r){const i=Aa;Aa=u;try{return r()}finally{Aa=i}}};return u}}let Aa=null;function Ze(s,n){if(Zs){let a=Zs.provides;const e=Zs.parent&&Zs.parent.provides;e===a&&(a=Zs.provides=Object.create(e)),a[s]=n}}function jn(s,n,a=!1){const e=ea();if(e||Aa){let t=Aa?Aa._context.provides:e?e.parent==null||e.ce?e.vnode.appContext&&e.vnode.appContext.provides:e.parent.provides:void 0;if(t&&s in t)return t[s];if(arguments.length>1)return a&&js(n)?n.call(e&&e.proxy):n}}function vr(){return!!(ea()||Aa)}const wr={},Cr=()=>Object.create(wr),kr=s=>Object.getPrototypeOf(s)===wr;function Ih(s,n,a,e=!1){const t={},l=Cr();s.propsDefaults=Object.create(null),xr(s,n,t,l);for(const o in s.propsOptions[0])o in t||(t[o]=void 0);a?s.props=e?t:Qc(t):s.type.props?s.props=t:s.props=l,s.attrs=l}function Nh(s,n,a,e){const{props:t,attrs:l,vnode:{patchFlag:o}}=s,p=xs(t),[c]=s.propsOptions;let u=!1;if((e||o>0)&&!(o&16)){if(o&8){const r=s.vnode.dynamicProps;for(let i=0;i<r.length;i++){let d=r[i];if(Et(s.emitsOptions,d))continue;const h=n[d];if(c)if(Es(l,d))h!==l[d]&&(l[d]=h,u=!0);else{const y=Tn(d);t[y]=_l(c,p,y,h,s,!1)}else h!==l[d]&&(l[d]=h,u=!0)}}}else{xr(s,n,t,l)&&(u=!0);let r;for(const i in p)(!n||!Es(n,i)&&((r=ga(i))===i||!Es(n,r)))&&(c?a&&(a[i]!==void 0||a[r]!==void 0)&&(t[i]=_l(c,p,i,void 0,s,!0)):delete t[i]);if(l!==p)for(const i in l)(!n||!Es(n,i))&&(delete l[i],u=!0)}u&&Xn(s.attrs,"set","")}function xr(s,n,a,e){const[t,l]=s.propsOptions;let o=!1,p;if(n)for(let c in n){if(te(c))continue;const u=n[c];let r;t&&Es(t,r=Tn(c))?!l||!l.includes(r)?a[r]=u:(p||(p={}))[r]=u:Et(s.emitsOptions,c)||(!(c in e)||u!==e[c])&&(e[c]=u,o=!0)}if(l){const c=xs(a),u=p||Ss;for(let r=0;r<l.length;r++){const i=l[r];a[i]=_l(t,c,i,u[i],s,!Es(u,i))}}return o}function _l(s,n,a,e,t,l){const o=s[a];if(o!=null){const p=Es(o,"default");if(p&&e===void 0){const c=o.default;if(o.type!==Function&&!o.skipFactory&&js(c)){const{propsDefaults:u}=t;if(a in u)e=u[a];else{const r=De(t);e=u[a]=c.call(null,n),r()}}else e=c;t.ce&&t.ce._setProp(a,e)}o[0]&&(l&&!p?e=!1:o[1]&&(e===""||e===ga(a))&&(e=!0))}return e}const Fh=new WeakMap;function Pr(s,n,a=!1){const e=a?Fh:n.propsCache,t=e.get(s);if(t)return t;const l=s.props,o={},p=[];let c=!1;if(!js(s)){const r=i=>{c=!0;const[d,h]=Pr(i,n,!0);Ks(o,d),h&&p.push(...h)};!a&&n.mixins.length&&n.mixins.forEach(r),s.extends&&r(s.extends),s.mixins&&s.mixins.forEach(r)}if(!l&&!c)return Ds(s)&&e.set(s,Ba),Ba;if(ds(l))for(let r=0;r<l.length;r++){const i=Tn(l[r]);Ho(i)&&(o[i]=Ss)}else if(l)for(const r in l){const i=Tn(r);if(Ho(i)){const d=l[r],h=o[i]=ds(d)||js(d)?{type:d}:Ks({},d),y=h.type;let j=!1,T=!0;if(ds(y))for(let g=0;g<y.length;++g){const m=y[g],w=js(m)&&m.name;if(w==="Boolean"){j=!0;break}else w==="String"&&(T=!1)}else j=js(y)&&y.name==="Boolean";h[0]=j,h[1]=T,(j||Es(h,"default"))&&p.push(i)}}const u=[o,p];return Ds(s)&&e.set(s,u),u}function Ho(s){return s[0]!=="$"&&!te(s)}const po=s=>s==="_"||s==="_ctx"||s==="$stable",co=s=>ds(s)?s.map(Fn):[Fn(s)],$h=(s,n,a)=>{if(n._n)return n;const e=vn((...t)=>co(n(...t)),a);return e._c=!1,e},Er=(s,n,a)=>{const e=s._ctx;for(const t in s){if(po(t))continue;const l=s[t];if(js(l))n[t]=$h(t,l,e);else if(l!=null){const o=co(l);n[t]=()=>o}}},Sr=(s,n)=>{const a=co(n);s.slots.default=()=>a},Tr=(s,n,a)=>{for(const e in n)(a||!po(e))&&(s[e]=n[e])},Bh=(s,n,a)=>{const e=s.slots=Cr();if(s.vnode.shapeFlag&32){const t=n._;t?(Tr(e,n,a),a&&Rc(e,"_",t,!0)):Er(n,e)}else n&&Sr(s,n)},qh=(s,n,a)=>{const{vnode:e,slots:t}=s;let l=!0,o=Ss;if(e.shapeFlag&32){const p=n._;p?a&&p===1?l=!1:Tr(t,n,a):(l=!n.$stable,Er(n,t)),o=n}else n&&(Sr(s,n),o={default:1});if(l)for(const p in t)!po(p)&&o[p]==null&&delete t[p]},rn=ed;function zh(s){return Vh(s)}function Vh(s,n){const a=kt();a.__VUE__=!0;const{insert:e,remove:t,patchProp:l,createElement:o,createText:p,createComment:c,setText:u,setElementText:r,parentNode:i,nextSibling:d,setScopeId:h=$n,insertStaticContent:y}=s,j=(C,P,R,$=null,V=null,q=null,b=void 0,v=null,A=!!P.dynamicChildren)=>{if(C===P)return;C&&!Ea(C,P)&&($=O(C),zs(C,V,q,!0),C=null),P.patchFlag===-2&&(A=!1,P.dynamicChildren=null);const{type:M,ref:Z,shapeFlag:K}=P;switch(M){case Me:T(C,P,R,$);break;case dn:g(C,P,R,$);break;case st:C==null&&m(P,R,$,b);break;case bs:U(C,P,R,$,V,q,b,v,A);break;default:K&1?x(C,P,R,$,V,q,b,v,A):K&6?F(C,P,R,$,V,q,b,v,A):(K&64||K&128)&&M.process(C,P,R,$,V,q,b,v,A,Q)}Z!=null&&V?ce(Z,C&&C.ref,q,P||C,!P):Z==null&&C&&C.ref!=null&&ce(C.ref,null,q,C,!0)},T=(C,P,R,$)=>{if(C==null)e(P.el=p(P.children),R,$);else{const V=P.el=C.el;P.children!==C.children&&u(V,P.children)}},g=(C,P,R,$)=>{C==null?e(P.el=c(P.children||""),R,$):P.el=C.el},m=(C,P,R,$)=>{[C.el,C.anchor]=y(C.children,P,R,$,C.el,C.anchor)},w=({el:C,anchor:P},R,$)=>{let V;for(;C&&C!==P;)V=d(C),e(C,R,$),C=V;e(P,R,$)},f=({el:C,anchor:P})=>{let R;for(;C&&C!==P;)R=d(C),t(C),C=R;t(P)},x=(C,P,R,$,V,q,b,v,A)=>{if(P.type==="svg"?b="svg":P.type==="math"&&(b="mathml"),C==null)k(P,R,$,V,q,b,v,A);else{const M=C.el&&C.el._isVueCE?C.el:null;try{M&&M._beginPatch(),S(C,P,V,q,b,v,A)}finally{M&&M._endPatch()}}},k=(C,P,R,$,V,q,b,v)=>{let A,M;const{props:Z,shapeFlag:K,transition:ns,dirs:as}=C;if(A=C.el=o(C.type,q,Z&&Z.is,Z),K&8?r(A,C.children):K&16&&D(C.children,A,null,$,V,Wt(C,q),b,v),as&&wa(C,null,$,"created"),_(A,C,C.scopeId,b,$),Z){for(const z in Z)z!=="value"&&!te(z)&&l(A,z,null,Z[z],q,$);"value"in Z&&l(A,"value",null,Z.value,q),(M=Z.onVnodeBeforeMount)&&On(M,$,C)}as&&wa(C,null,$,"beforeMount");const L=Uh(V,ns);L&&ns.beforeEnter(A),e(A,P,R),((M=Z&&Z.onVnodeMounted)||L||as)&&rn(()=>{M&&On(M,$,C),L&&ns.enter(A),as&&wa(C,null,$,"mounted")},V)},_=(C,P,R,$,V)=>{if(R&&h(C,R),$)for(let q=0;q<$.length;q++)h(C,$[q]);if(V){let q=V.subTree;if(P===q||Mr(q.type)&&(q.ssContent===P||q.ssFallback===P)){const b=V.vnode;_(C,b,b.scopeId,b.slotScopeIds,V.parent)}}},D=(C,P,R,$,V,q,b,v,A=0)=>{for(let M=A;M<C.length;M++){const Z=C[M]=v?ra(C[M]):Fn(C[M]);j(null,Z,P,R,$,V,q,b,v)}},S=(C,P,R,$,V,q,b)=>{const v=P.el=C.el;let{patchFlag:A,dynamicChildren:M,dirs:Z}=P;A|=C.patchFlag&16;const K=C.props||Ss,ns=P.props||Ss;let as;if(R&&Ca(R,!1),(as=ns.onVnodeBeforeUpdate)&&On(as,R,P,C),Z&&wa(P,C,R,"beforeUpdate"),R&&Ca(R,!0),(K.innerHTML&&ns.innerHTML==null||K.textContent&&ns.textContent==null)&&r(v,""),M?I(C.dynamicChildren,M,v,R,$,Wt(P,V),q):b||us(C,P,v,null,R,$,Wt(P,V),q,!1),A>0){if(A&16)J(v,K,ns,R,V);else if(A&2&&K.class!==ns.class&&l(v,"class",null,ns.class,V),A&4&&l(v,"style",K.style,ns.style,V),A&8){const L=P.dynamicProps;for(let z=0;z<L.length;z++){const ls=L[z],vs=K[ls],Bs=ns[ls];(Bs!==vs||ls==="value")&&l(v,ls,vs,Bs,V,R)}}A&1&&C.children!==P.children&&r(v,P.children)}else!b&&M==null&&J(v,K,ns,R,V);((as=ns.onVnodeUpdated)||Z)&&rn(()=>{as&&On(as,R,P,C),Z&&wa(P,C,R,"updated")},$)},I=(C,P,R,$,V,q,b)=>{for(let v=0;v<P.length;v++){const A=C[v],M=P[v],Z=A.el&&(A.type===bs||!Ea(A,M)||A.shapeFlag&198)?i(A.el):R;j(A,M,Z,null,$,V,q,b,!0)}},J=(C,P,R,$,V)=>{if(P!==R){if(P!==Ss)for(const q in P)!te(q)&&!(q in R)&&l(C,q,P[q],null,V,$);for(const q in R){if(te(q))continue;const b=R[q],v=P[q];b!==v&&q!=="value"&&l(C,q,v,b,V,$)}"value"in R&&l(C,"value",P.value,R.value,V)}},U=(C,P,R,$,V,q,b,v,A)=>{const M=P.el=C?C.el:p(""),Z=P.anchor=C?C.anchor:p("");let{patchFlag:K,dynamicChildren:ns,slotScopeIds:as}=P;as&&(v=v?v.concat(as):as),C==null?(e(M,R,$),e(Z,R,$),D(P.children||[],R,Z,V,q,b,v,A)):K>0&&K&64&&ns&&C.dynamicChildren?(I(C.dynamicChildren,ns,R,V,q,b,v),(P.key!=null||V&&P===V.subTree)&&ro(C,P,!0)):us(C,P,R,Z,V,q,b,v,A)},F=(C,P,R,$,V,q,b,v,A)=>{P.slotScopeIds=v,C==null?P.shapeFlag&512?V.ctx.activate(P,R,$,b,A):G(P,R,$,V,q,b,A):is(C,P,A)},G=(C,P,R,$,V,q,b)=>{const v=C.component=id(C,$,V);if(Re(C)&&(v.ctx.renderer=Q),ud(v,!1,b),v.asyncDep){if(V&&V.registerDep(v,ps,b),!C.el){const A=v.subTree=ws(dn);g(null,A,P,R),C.placeholder=A.el}}else ps(v,C,P,R,V,q,b)},is=(C,P,R)=>{const $=P.component=C.component;if(nd(C,P,R))if($.asyncDep&&!$.asyncResolved){Y($,P,R);return}else $.next=P,$.update();else P.el=C.el,$.vnode=P},ps=(C,P,R,$,V,q,b)=>{const v=()=>{if(C.isMounted){let{next:K,bu:ns,u:as,parent:L,vnode:z}=C;{const Ys=Ar(C);if(Ys){K&&(K.el=z.el,Y(C,K,b)),Ys.asyncDep.then(()=>{C.isUnmounted||v()});return}}let ls=K,vs;Ca(C,!1),K?(K.el=z.el,Y(C,K,b)):K=z,ns&&Ye(ns),(vs=K.props&&K.props.onVnodeBeforeUpdate)&&On(vs,L,K,z),Ca(C,!0);const Bs=Wo(C),on=C.subTree;C.subTree=Bs,j(on,Bs,i(on.el),O(on),C,V,q),K.el=Bs.el,ls===null&&ad(C,Bs.el),as&&rn(as,V),(vs=K.props&&K.props.onVnodeUpdated)&&rn(()=>On(vs,L,K,z),V)}else{let K;const{el:ns,props:as}=P,{bm:L,m:z,parent:ls,root:vs,type:Bs}=C,on=re(P);Ca(C,!1),L&&Ye(L),!on&&(K=as&&as.onVnodeBeforeMount)&&On(K,ls,P),Ca(C,!0);{vs.ce&&vs.ce._def.shadowRoot!==!1&&vs.ce._injectChildStyle(Bs);const Ys=C.subTree=Wo(C);j(null,Ys,R,$,C,V,q),P.el=Ys.el}if(z&&rn(z,V),!on&&(K=as&&as.onVnodeMounted)){const Ys=P;rn(()=>On(K,ls,Ys),V)}(P.shapeFlag&256||ls&&re(ls.vnode)&&ls.vnode.shapeFlag&256)&&C.a&&rn(C.a,V),C.isMounted=!0,P=R=$=null}};C.scope.on();const A=C.effect=new Nc(v);C.scope.off();const M=C.update=A.run.bind(A),Z=C.job=A.runIfDirty.bind(A);Z.i=C,Z.id=C.uid,A.scheduler=()=>no(Z),Ca(C,!0),M()},Y=(C,P,R)=>{P.component=C;const $=C.vnode.props;C.vnode=P,C.next=null,Nh(C,P.props,$,R),qh(C,P.children,R),Zn(),Mo(C),sa()},us=(C,P,R,$,V,q,b,v,A=!1)=>{const M=C&&C.children,Z=C?C.shapeFlag:0,K=P.children,{patchFlag:ns,shapeFlag:as}=P;if(ns>0){if(ns&128){Xs(M,K,R,$,V,q,b,v,A);return}else if(ns&256){Ns(M,K,R,$,V,q,b,v,A);return}}as&8?(Z&16&&_s(M,V,q),K!==M&&r(R,K)):Z&16?as&16?Xs(M,K,R,$,V,q,b,v,A):_s(M,V,q,!0):(Z&8&&r(R,""),as&16&&D(K,R,$,V,q,b,v,A))},Ns=(C,P,R,$,V,q,b,v,A)=>{C=C||Ba,P=P||Ba;const M=C.length,Z=P.length,K=Math.min(M,Z);let ns;for(ns=0;ns<K;ns++){const as=P[ns]=A?ra(P[ns]):Fn(P[ns]);j(C[ns],as,R,null,V,q,b,v,A)}M>Z?_s(C,V,q,!0,!1,K):D(P,R,$,V,q,b,v,A,K)},Xs=(C,P,R,$,V,q,b,v,A)=>{let M=0;const Z=P.length;let K=C.length-1,ns=Z-1;for(;M<=K&&M<=ns;){const as=C[M],L=P[M]=A?ra(P[M]):Fn(P[M]);if(Ea(as,L))j(as,L,R,null,V,q,b,v,A);else break;M++}for(;M<=K&&M<=ns;){const as=C[K],L=P[ns]=A?ra(P[ns]):Fn(P[ns]);if(Ea(as,L))j(as,L,R,null,V,q,b,v,A);else break;K--,ns--}if(M>K){if(M<=ns){const as=ns+1,L=as<Z?P[as].el:$;for(;M<=ns;)j(null,P[M]=A?ra(P[M]):Fn(P[M]),R,L,V,q,b,v,A),M++}}else if(M>ns)for(;M<=K;)zs(C[M],V,q,!0),M++;else{const as=M,L=M,z=new Map;for(M=L;M<=ns;M++){const an=P[M]=A?ra(P[M]):Fn(P[M]);an.key!=null&&z.set(an.key,M)}let ls,vs=0;const Bs=ns-L+1;let on=!1,Ys=0;const fn=new Array(Bs);for(M=0;M<Bs;M++)fn[M]=0;for(M=as;M<=K;M++){const an=C[M];if(vs>=Bs){zs(an,V,q,!0);continue}let Dn;if(an.key!=null)Dn=z.get(an.key);else for(ls=L;ls<=ns;ls++)if(fn[ls-L]===0&&Ea(an,P[ls])){Dn=ls;break}Dn===void 0?zs(an,V,q,!0):(fn[Dn-L]=M+1,Dn>=Ys?Ys=Dn:on=!0,j(an,P[Dn],R,null,V,q,b,v,A),vs++)}const Fe=on?Hh(fn):Ba;for(ls=Fe.length-1,M=Bs-1;M>=0;M--){const an=L+M,Dn=P[an],Co=P[an+1],ko=an+1<Z?Co.el||Co.placeholder:$;fn[M]===0?j(null,Dn,R,ko,V,q,b,v,A):on&&(ls<0||M!==Fe[ls]?Rs(Dn,R,ko,2):ls--)}}},Rs=(C,P,R,$,V=null)=>{const{el:q,type:b,transition:v,children:A,shapeFlag:M}=C;if(M&6){Rs(C.component.subTree,P,R,$);return}if(M&128){C.suspense.move(P,R,$);return}if(M&64){b.move(C,P,R,Q);return}if(b===bs){e(q,P,R);for(let K=0;K<A.length;K++)Rs(A[K],P,R,$);e(C.anchor,P,R);return}if(b===st){w(C,P,R);return}if($!==2&&M&1&&v)if($===0)v.beforeEnter(q),e(q,P,R),rn(()=>v.enter(q),V);else{const{leave:K,delayLeave:ns,afterLeave:as}=v,L=()=>{C.ctx.isUnmounted?t(q):e(q,P,R)},z=()=>{q._isLeaving&&q[Kn](!0),K(q,()=>{L(),as&&as()})};ns?ns(q,L,z):z()}else e(q,P,R)},zs=(C,P,R,$=!1,V=!1)=>{const{type:q,props:b,ref:v,children:A,dynamicChildren:M,shapeFlag:Z,patchFlag:K,dirs:ns,cacheIndex:as}=C;if(K===-2&&(V=!1),v!=null&&(Zn(),ce(v,null,R,C,!0),sa()),as!=null&&(P.renderCache[as]=void 0),Z&256){P.ctx.deactivate(C);return}const L=Z&1&&ns,z=!re(C);let ls;if(z&&(ls=b&&b.onVnodeBeforeUnmount)&&On(ls,P,C),Z&6)fs(C.component,R,$);else{if(Z&128){C.suspense.unmount(R,$);return}L&&wa(C,null,P,"beforeUnmount"),Z&64?C.type.remove(C,P,R,Q,$):M&&!M.hasOnce&&(q!==bs||K>0&&K&64)?_s(M,P,R,!1,!0):(q===bs&&K&384||!V&&Z&16)&&_s(A,P,R),$&&ts(C)}(z&&(ls=b&&b.onVnodeUnmounted)||L)&&rn(()=>{ls&&On(ls,P,C),L&&wa(C,null,P,"unmounted")},R)},ts=C=>{const{type:P,el:R,anchor:$,transition:V}=C;if(P===bs){rs(R,$);return}if(P===st){f(C);return}const q=()=>{t(R),V&&!V.persisted&&V.afterLeave&&V.afterLeave()};if(C.shapeFlag&1&&V&&!V.persisted){const{leave:b,delayLeave:v}=V,A=()=>b(R,q);v?v(C.el,q,A):A()}else q()},rs=(C,P)=>{let R;for(;C!==P;)R=d(C),t(C),C=R;t(P)},fs=(C,P,R)=>{const{bum:$,scope:V,job:q,subTree:b,um:v,m:A,a:M}=C;Go(A),Go(M),$&&Ye($),V.stop(),q&&(q.flags|=8,zs(b,C,P,R)),v&&rn(v,P),rn(()=>{C.isUnmounted=!0},P)},_s=(C,P,R,$=!1,V=!1,q=0)=>{for(let b=q;b<C.length;b++)zs(C[b],P,R,$,V)},O=C=>{if(C.shapeFlag&6)return O(C.component.subTree);if(C.shapeFlag&128)return C.suspense.next();const P=d(C.anchor||C.el),R=P&&P[lr];return R?d(R):P};let W=!1;const H=(C,P,R)=>{C==null?P._vnode&&zs(P._vnode,null,null,!0):j(P._vnode||null,C,P,null,null,null,R),P._vnode=C,W||(W=!0,Mo(),ar(),W=!1)},Q={p:j,um:zs,m:Rs,r:ts,mt:G,mc:D,pc:us,pbc:I,n:O,o:s};return{render:H,hydrate:void 0,createApp:Oh(H)}}function Wt({type:s,props:n},a){return a==="svg"&&s==="foreignObject"||a==="mathml"&&s==="annotation-xml"&&n&&n.encoding&&n.encoding.includes("html")?void 0:a}function Ca({effect:s,job:n},a){a?(s.flags|=32,n.flags|=4):(s.flags&=-33,n.flags&=-5)}function Uh(s,n){return(!s||s&&!s.pendingBranch)&&n&&!n.persisted}function ro(s,n,a=!1){const e=s.children,t=n.children;if(ds(e)&&ds(t))for(let l=0;l<e.length;l++){const o=e[l];let p=t[l];p.shapeFlag&1&&!p.dynamicChildren&&((p.patchFlag<=0||p.patchFlag===32)&&(p=t[l]=ra(t[l]),p.el=o.el),!a&&p.patchFlag!==-2&&ro(o,p)),p.type===Me&&p.patchFlag!==-1&&(p.el=o.el),p.type===dn&&!p.el&&(p.el=o.el)}}function Hh(s){const n=s.slice(),a=[0];let e,t,l,o,p;const c=s.length;for(e=0;e<c;e++){const u=s[e];if(u!==0){if(t=a[a.length-1],s[t]<u){n[e]=t,a.push(e);continue}for(l=0,o=a.length-1;l<o;)p=l+o>>1,s[a[p]]<u?l=p+1:o=p;u<s[a[l]]&&(l>0&&(n[e]=a[l-1]),a[l]=e)}}for(l=a.length,o=a[l-1];l-- >0;)a[l]=o,o=n[o];return a}function Ar(s){const n=s.subTree.component;if(n)return n.asyncDep&&!n.asyncResolved?n:Ar(n)}function Go(s){if(s)for(let n=0;n<s.length;n++)s[n].flags|=8}const Gh=Symbol.for("v-scx"),Wh=()=>jn(Gh);function Kh(s,n){return io(s,null,n)}function Hs(s,n,a){return io(s,n,a)}function io(s,n,a=Ss){const{immediate:e,deep:t,flush:l,once:o}=a,p=Ks({},a),c=n&&e||!n&&l!=="post";let u;if(Ha){if(l==="sync"){const h=Wh();u=h.__watcherHandles||(h.__watcherHandles=[])}else if(!c){const h=()=>{};return h.stop=$n,h.resume=$n,h.pause=$n,h}}const r=Zs;p.call=(h,y,j)=>Ln(h,r,y,j);let i=!1;l==="post"?p.scheduler=h=>{rn(h,r&&r.suspense)}:l!=="sync"&&(i=!0,p.scheduler=(h,y)=>{y?h():no(h)}),p.augmentJob=h=>{n&&(h.flags|=4),i&&(h.flags|=2,r&&(h.id=r.uid,h.i=r))};const d=ch(s,n,p);return Ha&&(u?u.push(d):c&&d()),d}function Xh(s,n,a){const e=this.proxy,t=$s(s)?s.includes(".")?Rr(e,s):()=>e[s]:s.bind(e,e);let l;js(n)?l=n:(l=n.handler,a=n);const o=De(this),p=io(t,l.bind(e),a);return o(),p}function Rr(s,n){const a=n.split(".");return()=>{let e=s;for(let t=0;t<a.length&&e;t++)e=e[a[t]];return e}}const Yh=(s,n)=>n==="modelValue"||n==="model-value"?s.modelModifiers:s[`${n}Modifiers`]||s[`${Tn(n)}Modifiers`]||s[`${ga(n)}Modifiers`];function Qh(s,n,...a){if(s.isUnmounted)return;const e=s.vnode.props||Ss;let t=a;const l=n.startsWith("update:"),o=l&&Yh(e,n.slice(7));o&&(o.trim&&(t=a.map(r=>$s(r)?r.trim():r)),o.number&&(t=a.map(Ul)));let p,c=e[p=Bt(n)]||e[p=Bt(Tn(n))];!c&&l&&(c=e[p=Bt(ga(n))]),c&&Ln(c,s,6,t);const u=e[p+"Once"];if(u){if(!s.emitted)s.emitted={};else if(s.emitted[p])return;s.emitted[p]=!0,Ln(u,s,6,t)}}const Jh=new WeakMap;function Lr(s,n,a=!1){const e=a?Jh:n.emitsCache,t=e.get(s);if(t!==void 0)return t;const l=s.emits;let o={},p=!1;if(!js(s)){const c=u=>{const r=Lr(u,n,!0);r&&(p=!0,Ks(o,r))};!a&&n.mixins.length&&n.mixins.forEach(c),s.extends&&c(s.extends),s.mixins&&s.mixins.forEach(c)}return!l&&!p?(Ds(s)&&e.set(s,null),null):(ds(l)?l.forEach(c=>o[c]=null):Ks(o,l),Ds(s)&&e.set(s,o),o)}function Et(s,n){return!s||!yt(n)?!1:(n=n.slice(2).replace(/Once$/,""),Es(s,n[0].toLowerCase()+n.slice(1))||Es(s,ga(n))||Es(s,n))}function Wo(s){const{type:n,vnode:a,proxy:e,withProxy:t,propsOptions:[l],slots:o,attrs:p,emit:c,render:u,renderCache:r,props:i,data:d,setupState:h,ctx:y,inheritAttrs:j}=s,T=rt(s);let g,m;try{if(a.shapeFlag&4){const f=t||e,x=f;g=Fn(u.call(x,f,r,i,h,d,y)),m=p}else{const f=n;g=Fn(f.length>1?f(i,{attrs:p,slots:o,emit:c}):f(i,null)),m=n.props?p:Zh(p)}}catch(f){ue.length=0,Ae(f,s,1),g=ws(dn)}let w=g;if(m&&j!==!1){const f=Object.keys(m),{shapeFlag:x}=w;f.length&&x&7&&(l&&f.some(ql)&&(m=sd(m,l)),w=ja(w,m,!1,!0))}return a.dirs&&(w=ja(w,null,!1,!0),w.dirs=w.dirs?w.dirs.concat(a.dirs):a.dirs),a.transition&&_e(w,a.transition),g=w,rt(T),g}const Zh=s=>{let n;for(const a in s)(a==="class"||a==="style"||yt(a))&&((n||(n={}))[a]=s[a]);return n},sd=(s,n)=>{const a={};for(const e in s)(!ql(e)||!(e.slice(9)in n))&&(a[e]=s[e]);return a};function nd(s,n,a){const{props:e,children:t,component:l}=s,{props:o,children:p,patchFlag:c}=n,u=l.emitsOptions;if(n.dirs||n.transition)return!0;if(a&&c>=0){if(c&1024)return!0;if(c&16)return e?Ko(e,o,u):!!o;if(c&8){const r=n.dynamicProps;for(let i=0;i<r.length;i++){const d=r[i];if(o[d]!==e[d]&&!Et(u,d))return!0}}}else return(t||p)&&(!p||!p.$stable)?!0:e===o?!1:e?o?Ko(e,o,u):!0:!!o;return!1}function Ko(s,n,a){const e=Object.keys(n);if(e.length!==Object.keys(s).length)return!0;for(let t=0;t<e.length;t++){const l=e[t];if(n[l]!==s[l]&&!Et(a,l))return!0}return!1}function ad({vnode:s,parent:n},a){for(;n;){const e=n.subTree;if(e.suspense&&e.suspense.activeBranch===s&&(e.el=s.el),e===s)(s=n.vnode).el=a,n=n.parent;else break}}const Mr=s=>s.__isSuspense;function ed(s,n){n&&n.pendingBranch?ds(s)?n.effects.push(...s):n.effects.push(s):uh(s)}const bs=Symbol.for("v-fgt"),Me=Symbol.for("v-txt"),dn=Symbol.for("v-cmt"),st=Symbol.for("v-stc"),ue=[];let bn=null;function N(s=!1){ue.push(bn=s?null:[])}function td(){ue.pop(),bn=ue[ue.length-1]||null}let ye=1;function ht(s,n=!1){ye+=s,s<0&&bn&&n&&(bn.hasOnce=!0)}function Dr(s){return s.dynamicChildren=ye>0?bn||Ba:null,td(),ye>0&&bn&&bn.push(s),s}function B(s,n,a,e,t,l){return Dr(E(s,n,a,e,t,l,!0))}function Sa(s,n,a,e,t){return Dr(ws(s,n,a,e,t,!0))}function dt(s){return s?s.__v_isVNode===!0:!1}function Ea(s,n){return s.type===n.type&&s.key===n.key}const Or=({key:s})=>s??null,nt=({ref:s,ref_key:n,ref_for:a})=>(typeof s=="number"&&(s=""+s),s!=null?$s(s)||Fs(s)||js(s)?{i:gn,r:s,k:n,f:!!a}:s:null);function E(s,n=null,a=null,e=0,t=null,l=s===bs?0:1,o=!1,p=!1){const c={__v_isVNode:!0,__v_skip:!0,type:s,props:n,key:n&&Or(n),ref:n&&nt(n),scopeId:tr,slotScopeIds:null,children:a,component:null,suspense:null,ssContent:null,ssFallback:null,dirs:null,transition:null,el:null,anchor:null,target:null,targetStart:null,targetAnchor:null,staticCount:0,shapeFlag:l,patchFlag:e,dynamicProps:t,dynamicChildren:null,appContext:null,ctx:gn};return p?(uo(c,a),l&128&&s.normalize(c)):a&&(c.shapeFlag|=$s(a)?8:16),ye>0&&!o&&bn&&(c.patchFlag>0||l&6)&&c.patchFlag!==32&&bn.push(c),c}const ws=ld;function ld(s,n=null,a=null,e=0,t=null,l=!1){if((!s||s===gr)&&(s=dn),dt(s)){const p=ja(s,n,!0);return a&&uo(p,a),ye>0&&!l&&bn&&(p.shapeFlag&6?bn[bn.indexOf(s)]=p:bn.push(p)),p.patchFlag=-2,p}if(fd(s)&&(s=s.__vccOpts),n){n=od(n);let{class:p,style:c}=n;p&&!$s(p)&&(n.class=Us(p)),Ds(c)&&(Jl(c)&&!ds(c)&&(c=Ks({},c)),n.style=ba(c))}const o=$s(s)?1:Mr(s)?128:or(s)?64:Ds(s)?4:js(s)?2:0;return E(s,n,a,e,t,o,l,!0)}function od(s){return s?Jl(s)||kr(s)?Ks({},s):s:null}function ja(s,n,a=!1,e=!1){const{props:t,ref:l,patchFlag:o,children:p,transition:c}=s,u=n?pd(t||{},n):t,r={__v_isVNode:!0,__v_skip:!0,type:s.type,props:u,key:u&&Or(u),ref:n&&n.ref?a&&l?ds(l)?l.concat(nt(n)):[l,nt(n)]:nt(n):l,scopeId:s.scopeId,slotScopeIds:s.slotScopeIds,children:p,target:s.target,targetStart:s.targetStart,targetAnchor:s.targetAnchor,staticCount:s.staticCount,shapeFlag:s.shapeFlag,patchFlag:n&&s.type!==bs?o===-1?16:o|16:o,dynamicProps:s.dynamicProps,dynamicChildren:s.dynamicChildren,appContext:s.appContext,dirs:s.dirs,transition:c,component:s.component,suspense:s.suspense,ssContent:s.ssContent&&ja(s.ssContent),ssFallback:s.ssFallback&&ja(s.ssFallback),placeholder:s.placeholder,el:s.el,anchor:s.anchor,ctx:s.ctx,ce:s.ce};return c&&e&&_e(r,c.clone(r)),r}function As(s=" ",n=0){return ws(Me,null,s,n)}function He(s,n){const a=ws(st,null,s);return a.staticCount=n,a}function hs(s="",n=!1){return n?(N(),Sa(dn,null,s)):ws(dn,null,s)}function Fn(s){return s==null||typeof s=="boolean"?ws(dn):ds(s)?ws(bs,null,s.slice()):dt(s)?ra(s):ws(Me,null,String(s))}function ra(s){return s.el===null&&s.patchFlag!==-1||s.memo?s:ja(s)}function uo(s,n){let a=0;const{shapeFlag:e}=s;if(n==null)n=null;else if(ds(n))a=16;else if(typeof n=="object")if(e&65){const t=n.default;t&&(t._c&&(t._d=!1),uo(s,t()),t._c&&(t._d=!0));return}else{a=32;const t=n._;!t&&!kr(n)?n._ctx=gn:t===3&&gn&&(gn.slots._===1?n._=1:(n._=2,s.patchFlag|=1024))}else js(n)?(n={default:n,_ctx:gn},a=32):(n=String(n),e&64?(a=16,n=[As(n)]):a=8);s.children=n,s.shapeFlag|=a}function pd(...s){const n={};for(let a=0;a<s.length;a++){const e=s[a];for(const t in e)if(t==="class")n.class!==e.class&&(n.class=Us([n.class,e.class]));else if(t==="style")n.style=ba([n.style,e.style]);else if(yt(t)){const l=n[t],o=e[t];o&&l!==o&&!(ds(l)&&l.includes(o))&&(n[t]=l?[].concat(l,o):o)}else t!==""&&(n[t]=e[t])}return n}function On(s,n,a,e=null){Ln(s,n,7,[a,e])}const cd=yr();let rd=0;function id(s,n,a){const e=s.type,t=(n?n.appContext:s.appContext)||cd,l={uid:rd++,vnode:s,type:e,parent:n,appContext:t,root:null,next:null,subTree:null,effect:null,update:null,job:null,scope:new Oc(!0),render:null,proxy:null,exposed:null,exposeProxy:null,withProxy:null,provides:n?n.provides:Object.create(t.provides),ids:n?n.ids:["",0,0],accessCache:null,renderCache:[],components:null,directives:null,propsOptions:Pr(e,t),emitsOptions:Lr(e,t),emit:null,emitted:null,propsDefaults:Ss,inheritAttrs:e.inheritAttrs,ctx:Ss,data:Ss,props:Ss,attrs:Ss,slots:Ss,refs:Ss,setupState:Ss,setupContext:null,suspense:a,suspenseId:a?a.pendingId:0,asyncDep:null,asyncResolved:!1,isMounted:!1,isUnmounted:!1,isDeactivated:!1,bc:null,c:null,bm:null,m:null,bu:null,u:null,um:null,bum:null,da:null,a:null,rtg:null,rtc:null,ec:null,sp:null};return l.ctx={_:l},l.root=n?n.root:l,l.emit=Qh.bind(null,l),s.ce&&s.ce(l),l}let Zs=null;const ea=()=>Zs||gn;let mt,yl;{const s=kt(),n=(a,e)=>{let t;return(t=s[a])||(t=s[a]=[]),t.push(e),l=>{t.length>1?t.forEach(o=>o(l)):t[0](l)}};mt=n("__VUE_INSTANCE_SETTERS__",a=>Zs=a),yl=n("__VUE_SSR_SETTERS__",a=>Ha=a)}const De=s=>{const n=Zs;return mt(s),s.scope.on(),()=>{s.scope.off(),mt(n)}},Xo=()=>{Zs&&Zs.scope.off(),mt(null)};function Ir(s){return s.vnode.shapeFlag&4}let Ha=!1;function ud(s,n=!1,a=!1){n&&yl(n);const{props:e,children:t}=s.vnode,l=Ir(s);Ih(s,e,l,n),Bh(s,t,a||n);const o=l?hd(s,n):void 0;return n&&yl(!1),o}function hd(s,n){const a=s.type;s.accessCache=Object.create(null),s.proxy=new Proxy(s.ctx,Sh);const{setup:e}=a;if(e){Zn();const t=s.setupContext=e.length>1?md(s):null,l=De(s),o=Te(e,s,0,[s.props,t]),p=Sc(o);if(sa(),l(),(p||s.sp)&&!re(s)&&ao(s),p){if(o.then(Xo,Xo),n)return o.then(c=>{Yo(s,c)}).catch(c=>{Ae(c,s,0)});s.asyncDep=o}else Yo(s,o)}else Nr(s)}function Yo(s,n,a){js(n)?s.type.__ssrInlineRender?s.ssrRender=n:s.render=n:Ds(n)&&(s.setupState=Zc(n)),Nr(s)}function Nr(s,n,a){const e=s.type;s.render||(s.render=e.render||$n);{const t=De(s);Zn();try{Th(s)}finally{sa(),t()}}}const dd={get(s,n){return ln(s,"get",""),s[n]}};function md(s){const n=a=>{s.exposed=a||{}};return{attrs:new Proxy(s.attrs,dd),slots:s.slots,emit:s.emit,expose:n}}function St(s){return s.exposed?s.exposeProxy||(s.exposeProxy=new Proxy(Zc(Zl(s.exposed)),{get(n,a){if(a in n)return n[a];if(a in ie)return ie[a](s)},has(n,a){return a in n||a in ie}})):s.proxy}function jd(s,n=!0){return js(s)?s.displayName||s.name:s.name||n&&s.__name}function fd(s){return js(s)&&"__vccOpts"in s}const ss=(s,n)=>oh(s,n,Ha);function Oe(s,n,a){try{ht(-1);const e=arguments.length;return e===2?Ds(n)&&!ds(n)?dt(n)?ws(s,null,[n]):ws(s,n):ws(s,null,n):(e>3?a=Array.prototype.slice.call(arguments,2):e===3&&dt(a)&&(a=[a]),ws(s,n,a))}finally{ht(1)}}const gd="3.5.24";/**
* @vue/runtime-dom v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/let vl;const Qo=typeof window<"u"&&window.trustedTypes;if(Qo)try{vl=Qo.createPolicy("vue",{createHTML:s=>s})}catch{}const Fr=vl?s=>vl.createHTML(s):s=>s,bd="http://www.w3.org/2000/svg",_d="http://www.w3.org/1998/Math/MathML",Wn=typeof document<"u"?document:null,Jo=Wn&&Wn.createElement("template"),yd={insert:(s,n,a)=>{n.insertBefore(s,a||null)},remove:s=>{const n=s.parentNode;n&&n.removeChild(s)},createElement:(s,n,a,e)=>{const t=n==="svg"?Wn.createElementNS(bd,s):n==="mathml"?Wn.createElementNS(_d,s):a?Wn.createElement(s,{is:a}):Wn.createElement(s);return s==="select"&&e&&e.multiple!=null&&t.setAttribute("multiple",e.multiple),t},createText:s=>Wn.createTextNode(s),createComment:s=>Wn.createComment(s),setText:(s,n)=>{s.nodeValue=n},setElementText:(s,n)=>{s.textContent=n},parentNode:s=>s.parentNode,nextSibling:s=>s.nextSibling,querySelector:s=>Wn.querySelector(s),setScopeId(s,n){s.setAttribute(n,"")},insertStaticContent(s,n,a,e,t,l){const o=a?a.previousSibling:n.lastChild;if(t&&(t===l||t.nextSibling))for(;n.insertBefore(t.cloneNode(!0),a),!(t===l||!(t=t.nextSibling)););else{Jo.innerHTML=Fr(e==="svg"?`<svg>${s}</svg>`:e==="mathml"?`<math>${s}</math>`:s);const p=Jo.content;if(e==="svg"||e==="mathml"){const c=p.firstChild;for(;c.firstChild;)p.appendChild(c.firstChild);p.removeChild(c)}n.insertBefore(p,a)}return[o?o.nextSibling:n.firstChild,a?a.previousSibling:n.lastChild]}},ta="transition",Ja="animation",ve=Symbol("_vtc"),$r={name:String,type:String,css:{type:Boolean,default:!0},duration:[String,Number,Object],enterFromClass:String,enterActiveClass:String,enterToClass:String,appearFromClass:String,appearActiveClass:String,appearToClass:String,leaveFromClass:String,leaveActiveClass:String,leaveToClass:String},vd=Ks({},rr,$r),wd=s=>(s.displayName="Transition",s.props=vd,s),Cd=wd((s,{slots:n})=>Oe(jh,kd(s),n)),ka=(s,n=[])=>{ds(s)?s.forEach(a=>a(...n)):s&&s(...n)},Zo=s=>s?ds(s)?s.some(n=>n.length>1):s.length>1:!1;function kd(s){const n={};for(const U in s)U in $r||(n[U]=s[U]);if(s.css===!1)return n;const{name:a="v",type:e,duration:t,enterFromClass:l=`${a}-enter-from`,enterActiveClass:o=`${a}-enter-active`,enterToClass:p=`${a}-enter-to`,appearFromClass:c=l,appearActiveClass:u=o,appearToClass:r=p,leaveFromClass:i=`${a}-leave-from`,leaveActiveClass:d=`${a}-leave-active`,leaveToClass:h=`${a}-leave-to`}=s,y=xd(t),j=y&&y[0],T=y&&y[1],{onBeforeEnter:g,onEnter:m,onEnterCancelled:w,onLeave:f,onLeaveCancelled:x,onBeforeAppear:k=g,onAppear:_=m,onAppearCancelled:D=w}=n,S=(U,F,G,is)=>{U._enterCancelled=is,xa(U,F?r:p),xa(U,F?u:o),G&&G()},I=(U,F)=>{U._isLeaving=!1,xa(U,i),xa(U,h),xa(U,d),F&&F()},J=U=>(F,G)=>{const is=U?_:m,ps=()=>S(F,U,G);ka(is,[F,ps]),sp(()=>{xa(F,U?c:l),Un(F,U?r:p),Zo(is)||np(F,e,j,ps)})};return Ks(n,{onBeforeEnter(U){ka(g,[U]),Un(U,l),Un(U,o)},onBeforeAppear(U){ka(k,[U]),Un(U,c),Un(U,u)},onEnter:J(!1),onAppear:J(!0),onLeave(U,F){U._isLeaving=!0;const G=()=>I(U,F);Un(U,i),U._enterCancelled?(Un(U,d),tp(U)):(tp(U),Un(U,d)),sp(()=>{U._isLeaving&&(xa(U,i),Un(U,h),Zo(f)||np(U,e,T,G))}),ka(f,[U,G])},onEnterCancelled(U){S(U,!1,void 0,!0),ka(w,[U])},onAppearCancelled(U){S(U,!0,void 0,!0),ka(D,[U])},onLeaveCancelled(U){I(U),ka(x,[U])}})}function xd(s){if(s==null)return null;if(Ds(s))return[Kt(s.enter),Kt(s.leave)];{const n=Kt(s);return[n,n]}}function Kt(s){return Pu(s)}function Un(s,n){n.split(/\s+/).forEach(a=>a&&s.classList.add(a)),(s[ve]||(s[ve]=new Set)).add(n)}function xa(s,n){n.split(/\s+/).forEach(e=>e&&s.classList.remove(e));const a=s[ve];a&&(a.delete(n),a.size||(s[ve]=void 0))}function sp(s){requestAnimationFrame(()=>{requestAnimationFrame(s)})}let Pd=0;function np(s,n,a,e){const t=s._endId=++Pd,l=()=>{t===s._endId&&e()};if(a!=null)return setTimeout(l,a);const{type:o,timeout:p,propCount:c}=Ed(s,n);if(!o)return e();const u=o+"end";let r=0;const i=()=>{s.removeEventListener(u,d),l()},d=h=>{h.target===s&&++r>=c&&i()};setTimeout(()=>{r<c&&i()},p+1),s.addEventListener(u,d)}function Ed(s,n){const a=window.getComputedStyle(s),e=y=>(a[y]||"").split(", "),t=e(`${ta}Delay`),l=e(`${ta}Duration`),o=ap(t,l),p=e(`${Ja}Delay`),c=e(`${Ja}Duration`),u=ap(p,c);let r=null,i=0,d=0;n===ta?o>0&&(r=ta,i=o,d=l.length):n===Ja?u>0&&(r=Ja,i=u,d=c.length):(i=Math.max(o,u),r=i>0?o>u?ta:Ja:null,d=r?r===ta?l.length:c.length:0);const h=r===ta&&/\b(?:transform|all)(?:,|$)/.test(e(`${ta}Property`).toString());return{type:r,timeout:i,propCount:d,hasTransform:h}}function ap(s,n){for(;s.length<n.length;)s=s.concat(s);return Math.max(...n.map((a,e)=>ep(a)+ep(s[e])))}function ep(s){return s==="auto"?0:Number(s.slice(0,-1).replace(",","."))*1e3}function tp(s){return(s?s.ownerDocument:document).body.offsetHeight}function Sd(s,n,a){const e=s[ve];e&&(n=(n?[n,...e]:[...e]).join(" ")),n==null?s.removeAttribute("class"):a?s.setAttribute("class",n):s.className=n}const jt=Symbol("_vod"),Br=Symbol("_vsh"),qr={name:"show",beforeMount(s,{value:n},{transition:a}){s[jt]=s.style.display==="none"?"":s.style.display,a&&n?a.beforeEnter(s):Za(s,n)},mounted(s,{value:n},{transition:a}){a&&n&&a.enter(s)},updated(s,{value:n,oldValue:a},{transition:e}){!n!=!a&&(e?n?(e.beforeEnter(s),Za(s,!0),e.enter(s)):e.leave(s,()=>{Za(s,!1)}):Za(s,n))},beforeUnmount(s,{value:n}){Za(s,n)}};function Za(s,n){s.style.display=n?s[jt]:"none",s[Br]=!n}const Td=Symbol(""),Ad=/(?:^|;)\s*display\s*:/;function Rd(s,n,a){const e=s.style,t=$s(a);let l=!1;if(a&&!t){if(n)if($s(n))for(const o of n.split(";")){const p=o.slice(0,o.indexOf(":")).trim();a[p]==null&&at(e,p,"")}else for(const o in n)a[o]==null&&at(e,o,"");for(const o in a)o==="display"&&(l=!0),at(e,o,a[o])}else if(t){if(n!==a){const o=e[Td];o&&(a+=";"+o),e.cssText=a,l=Ad.test(a)}}else n&&s.removeAttribute("style");jt in s&&(s[jt]=l?e.display:"",s[Br]&&(e.display="none"))}const lp=/\s*!important$/;function at(s,n,a){if(ds(a))a.forEach(e=>at(s,n,e));else if(a==null&&(a=""),n.startsWith("--"))s.setProperty(n,a);else{const e=Ld(s,n);lp.test(a)?s.setProperty(ga(e),a.replace(lp,""),"important"):s[e]=a}}const op=["Webkit","Moz","ms"],Xt={};function Ld(s,n){const a=Xt[n];if(a)return a;let e=Tn(n);if(e!=="filter"&&e in s)return Xt[n]=e;e=Ct(e);for(let t=0;t<op.length;t++){const l=op[t]+e;if(l in s)return Xt[n]=l}return n}const pp="http://www.w3.org/1999/xlink";function cp(s,n,a,e,t,l=Lu(n)){e&&n.startsWith("xlink:")?a==null?s.removeAttributeNS(pp,n.slice(6,n.length)):s.setAttributeNS(pp,n,a):a==null||l&&!Lc(a)?s.removeAttribute(n):s.setAttribute(n,l?"":fa(a)?String(a):a)}function rp(s,n,a,e,t){if(n==="innerHTML"||n==="textContent"){a!=null&&(s[n]=n==="innerHTML"?Fr(a):a);return}const l=s.tagName;if(n==="value"&&l!=="PROGRESS"&&!l.includes("-")){const p=l==="OPTION"?s.getAttribute("value")||"":s.value,c=a==null?s.type==="checkbox"?"on":"":String(a);(p!==c||!("_value"in s))&&(s.value=c),a==null&&s.removeAttribute(n),s._value=a;return}let o=!1;if(a===""||a==null){const p=typeof s[n];p==="boolean"?a=Lc(a):a==null&&p==="string"?(a="",o=!0):p==="number"&&(a=0,o=!0)}try{s[n]=a}catch{}o&&s.removeAttribute(t||n)}function Na(s,n,a,e){s.addEventListener(n,a,e)}function Md(s,n,a,e){s.removeEventListener(n,a,e)}const ip=Symbol("_vei");function Dd(s,n,a,e,t=null){const l=s[ip]||(s[ip]={}),o=l[n];if(e&&o)o.value=e;else{const[p,c]=Od(n);if(e){const u=l[n]=Fd(e,t);Na(s,p,u,c)}else o&&(Md(s,p,o,c),l[n]=void 0)}}const up=/(?:Once|Passive|Capture)$/;function Od(s){let n;if(up.test(s)){n={};let e;for(;e=s.match(up);)s=s.slice(0,s.length-e[0].length),n[e[0].toLowerCase()]=!0}return[s[2]===":"?s.slice(3):ga(s.slice(2)),n]}let Yt=0;const Id=Promise.resolve(),Nd=()=>Yt||(Id.then(()=>Yt=0),Yt=Date.now());function Fd(s,n){const a=e=>{if(!e._vts)e._vts=Date.now();else if(e._vts<=a.attached)return;Ln($d(e,a.value),n,5,[e])};return a.value=s,a.attached=Nd(),a}function $d(s,n){if(ds(n)){const a=s.stopImmediatePropagation;return s.stopImmediatePropagation=()=>{a.call(s),s._stopped=!0},n.map(e=>t=>!t._stopped&&e&&e(t))}else return n}const hp=s=>s.charCodeAt(0)===111&&s.charCodeAt(1)===110&&s.charCodeAt(2)>96&&s.charCodeAt(2)<123,Bd=(s,n,a,e,t,l)=>{const o=t==="svg";n==="class"?Sd(s,e,o):n==="style"?Rd(s,a,e):yt(n)?ql(n)||Dd(s,n,a,e,l):(n[0]==="."?(n=n.slice(1),!0):n[0]==="^"?(n=n.slice(1),!1):qd(s,n,e,o))?(rp(s,n,e),!s.tagName.includes("-")&&(n==="value"||n==="checked"||n==="selected")&&cp(s,n,e,o,l,n!=="value")):s._isVueCE&&(/[A-Z]/.test(n)||!$s(e))?rp(s,Tn(n),e,l,n):(n==="true-value"?s._trueValue=e:n==="false-value"&&(s._falseValue=e),cp(s,n,e,o))};function qd(s,n,a,e){if(e)return!!(n==="innerHTML"||n==="textContent"||n in s&&hp(n)&&js(a));if(n==="spellcheck"||n==="draggable"||n==="translate"||n==="autocorrect"||n==="sandbox"&&s.tagName==="IFRAME"||n==="form"||n==="list"&&s.tagName==="INPUT"||n==="type"&&s.tagName==="TEXTAREA")return!1;if(n==="width"||n==="height"){const t=s.tagName;if(t==="IMG"||t==="VIDEO"||t==="CANVAS"||t==="SOURCE")return!1}return hp(n)&&$s(a)?!1:n in s}const dp=s=>{const n=s.props["onUpdate:modelValue"]||!1;return ds(n)?a=>Ye(n,a):n};function zd(s){s.target.composing=!0}function mp(s){const n=s.target;n.composing&&(n.composing=!1,n.dispatchEvent(new Event("input")))}const Qt=Symbol("_assign");function jp(s,n,a){return n&&(s=s.trim()),a&&(s=Ul(s)),s}const kw={created(s,{modifiers:{lazy:n,trim:a,number:e}},t){s[Qt]=dp(t);const l=e||t.props&&t.props.type==="number";Na(s,n?"change":"input",o=>{o.target.composing||s[Qt](jp(s.value,a,l))}),(a||l)&&Na(s,"change",()=>{s.value=jp(s.value,a,l)}),n||(Na(s,"compositionstart",zd),Na(s,"compositionend",mp),Na(s,"change",mp))},mounted(s,{value:n}){s.value=n??""},beforeUpdate(s,{value:n,oldValue:a,modifiers:{lazy:e,trim:t,number:l}},o){if(s[Qt]=dp(o),s.composing)return;const p=(l||s.type==="number")&&!/^0\d/.test(s.value)?Ul(s.value):s.value,c=n??"";p!==c&&(document.activeElement===s&&s.type!=="range"&&(e&&n===a||t&&s.value.trim()===c)||(s.value=c))}},Vd=["ctrl","shift","alt","meta"],Ud={stop:s=>s.stopPropagation(),prevent:s=>s.preventDefault(),self:s=>s.target!==s.currentTarget,ctrl:s=>!s.ctrlKey,shift:s=>!s.shiftKey,alt:s=>!s.altKey,meta:s=>!s.metaKey,left:s=>"button"in s&&s.button!==0,middle:s=>"button"in s&&s.button!==1,right:s=>"button"in s&&s.button!==2,exact:(s,n)=>Vd.some(a=>s[`${a}Key`]&&!n.includes(a))},_n=(s,n)=>{const a=s._withMods||(s._withMods={}),e=n.join(".");return a[e]||(a[e]=((t,...l)=>{for(let o=0;o<n.length;o++){const p=Ud[n[o]];if(p&&p(t,n))return}return s(t,...l)}))},Hd={esc:"escape",space:" ",up:"arrow-up",left:"arrow-left",right:"arrow-right",down:"arrow-down",delete:"backspace"},fp=(s,n)=>{const a=s._withKeys||(s._withKeys={}),e=n.join(".");return a[e]||(a[e]=(t=>{if(!("key"in t))return;const l=ga(t.key);if(n.some(o=>o===l||Hd[o]===l))return s(t)}))},Gd=Ks({patchProp:Bd},yd);let gp;function Wd(){return gp||(gp=zh(Gd))}const Kd=((...s)=>{const n=Wd().createApp(...s),{mount:a}=n;return n.mount=e=>{const t=Yd(e);if(!t)return;const l=n._component;!js(l)&&!l.render&&!l.template&&(l.template=t.innerHTML),t.nodeType===1&&(t.textContent="");const o=a(t,!1,Xd(t));return t instanceof Element&&(t.removeAttribute("v-cloak"),t.setAttribute("data-v-app","")),o},n});function Xd(s){if(s instanceof SVGElement)return"svg";if(typeof MathMLElement=="function"&&s instanceof MathMLElement)return"mathml"}function Yd(s){return $s(s)?document.querySelector(s):s}const Qd=(s,n)=>Fs(n)?sh(n):n,zr="usehead";function Jd(s){return{install(a){a.config.globalProperties.$unhead=s,a.config.globalProperties.$head=s,a.provide(zr,s)}}.install}function Zd(){if(vr()){const s=jn(zr);if(s)return s}throw new Error("useHead() was called without provide context, ensure you call it through the setup() function.")}function Tt(s,n={}){const a=n.head||Zd();return a.ssr?a.push(s||{},n):sm(a,s,n)}function sm(s,n,a={}){const e=os(!1);let t;return Kh(()=>{const o=e.value?{}:tt(n,Qd);t?t.patch(o):t=s.push(o,a)}),ea()&&(qn(()=>{t.dispose()}),jr(()=>{e.value=!0}),mr(()=>{e.value=!1})),t}const nm={created(){let s=!1;const n=ea();if(!n)return;const a=n.type;!a||!("head"in a)||(s=typeof a.head=="function"?()=>a.head.call(n.proxy):a.head,s&&Tt(s))}};function am(s={}){const n=yu({domOptions:{render:vu(()=>xc(n),a=>setTimeout(a,0))},...s});return n.install=Jd(n),n}/*!
  * vue-router v4.5.1
  * (c) 2025 Eduardo San Martin Morote
  * @license MIT
  */const Fa=typeof document<"u";function Vr(s){return typeof s=="object"||"displayName"in s||"props"in s||"__vccOpts"in s}function em(s){return s.__esModule||s[Symbol.toStringTag]==="Module"||s.default&&Vr(s.default)}const Ps=Object.assign;function Jt(s,n){const a={};for(const e in n){const t=n[e];a[e]=Mn(t)?t.map(s):s(t)}return a}const he=()=>{},Mn=Array.isArray,Ur=/#/g,tm=/&/g,lm=/\//g,om=/=/g,pm=/\?/g,Hr=/\+/g,cm=/%5B/g,rm=/%5D/g,Gr=/%5E/g,im=/%60/g,Wr=/%7B/g,um=/%7C/g,Kr=/%7D/g,hm=/%20/g;function ho(s){return encodeURI(""+s).replace(um,"|").replace(cm,"[").replace(rm,"]")}function dm(s){return ho(s).replace(Wr,"{").replace(Kr,"}").replace(Gr,"^")}function wl(s){return ho(s).replace(Hr,"%2B").replace(hm,"+").replace(Ur,"%23").replace(tm,"%26").replace(im,"`").replace(Wr,"{").replace(Kr,"}").replace(Gr,"^")}function mm(s){return wl(s).replace(om,"%3D")}function jm(s){return ho(s).replace(Ur,"%23").replace(pm,"%3F")}function fm(s){return s==null?"":jm(s).replace(lm,"%2F")}function we(s){try{return decodeURIComponent(""+s)}catch{}return""+s}const gm=/\/$/,bm=s=>s.replace(gm,"");function Zt(s,n,a="/"){let e,t={},l="",o="";const p=n.indexOf("#");let c=n.indexOf("?");return p<c&&p>=0&&(c=-1),c>-1&&(e=n.slice(0,c),l=n.slice(c+1,p>-1?p:n.length),t=s(l)),p>-1&&(e=e||n.slice(0,p),o=n.slice(p,n.length)),e=wm(e??n,a),{fullPath:e+(l&&"?")+l+o,path:e,query:t,hash:we(o)}}function _m(s,n){const a=n.query?s(n.query):"";return n.path+(a&&"?")+a+(n.hash||"")}function bp(s,n){return!n||!s.toLowerCase().startsWith(n.toLowerCase())?s:s.slice(n.length)||"/"}function ym(s,n,a){const e=n.matched.length-1,t=a.matched.length-1;return e>-1&&e===t&&Ga(n.matched[e],a.matched[t])&&Xr(n.params,a.params)&&s(n.query)===s(a.query)&&n.hash===a.hash}function Ga(s,n){return(s.aliasOf||s)===(n.aliasOf||n)}function Xr(s,n){if(Object.keys(s).length!==Object.keys(n).length)return!1;for(const a in s)if(!vm(s[a],n[a]))return!1;return!0}function vm(s,n){return Mn(s)?_p(s,n):Mn(n)?_p(n,s):s===n}function _p(s,n){return Mn(n)?s.length===n.length&&s.every((a,e)=>a===n[e]):s.length===1&&s[0]===n}function wm(s,n){if(s.startsWith("/"))return s;if(!s)return n;const a=n.split("/"),e=s.split("/"),t=e[e.length-1];(t===".."||t===".")&&e.push("");let l=a.length-1,o,p;for(o=0;o<e.length;o++)if(p=e[o],p!==".")if(p==="..")l>1&&l--;else break;return a.slice(0,l).join("/")+"/"+e.slice(o).join("/")}const la={path:"/",name:void 0,params:{},query:{},hash:"",fullPath:"/",matched:[],meta:{},redirectedFrom:void 0};var Ce;(function(s){s.pop="pop",s.push="push"})(Ce||(Ce={}));var de;(function(s){s.back="back",s.forward="forward",s.unknown=""})(de||(de={}));function Cm(s){if(!s)if(Fa){const n=document.querySelector("base");s=n&&n.getAttribute("href")||"/",s=s.replace(/^\w+:\/\/[^\/]+/,"")}else s="/";return s[0]!=="/"&&s[0]!=="#"&&(s="/"+s),bm(s)}const km=/^[^#]+#/;function xm(s,n){return s.replace(km,"#")+n}function Pm(s,n){const a=document.documentElement.getBoundingClientRect(),e=s.getBoundingClientRect();return{behavior:n.behavior,left:e.left-a.left-(n.left||0),top:e.top-a.top-(n.top||0)}}const At=()=>({left:window.scrollX,top:window.scrollY});function Em(s){let n;if("el"in s){const a=s.el,e=typeof a=="string"&&a.startsWith("#"),t=typeof a=="string"?e?document.getElementById(a.slice(1)):document.querySelector(a):a;if(!t)return;n=Pm(t,s)}else n=s;"scrollBehavior"in document.documentElement.style?window.scrollTo(n):window.scrollTo(n.left!=null?n.left:window.scrollX,n.top!=null?n.top:window.scrollY)}function yp(s,n){return(history.state?history.state.position-n:-1)+s}const Cl=new Map;function Sm(s,n){Cl.set(s,n)}function Tm(s){const n=Cl.get(s);return Cl.delete(s),n}let Am=()=>location.protocol+"//"+location.host;function Yr(s,n){const{pathname:a,search:e,hash:t}=n,l=s.indexOf("#");if(l>-1){let p=t.includes(s.slice(l))?s.slice(l).length:1,c=t.slice(p);return c[0]!=="/"&&(c="/"+c),bp(c,"")}return bp(a,s)+e+t}function Rm(s,n,a,e){let t=[],l=[],o=null;const p=({state:d})=>{const h=Yr(s,location),y=a.value,j=n.value;let T=0;if(d){if(a.value=h,n.value=d,o&&o===y){o=null;return}T=j?d.position-j.position:0}else e(h);t.forEach(g=>{g(a.value,y,{delta:T,type:Ce.pop,direction:T?T>0?de.forward:de.back:de.unknown})})};function c(){o=a.value}function u(d){t.push(d);const h=()=>{const y=t.indexOf(d);y>-1&&t.splice(y,1)};return l.push(h),h}function r(){const{history:d}=window;d.state&&d.replaceState(Ps({},d.state,{scroll:At()}),"")}function i(){for(const d of l)d();l=[],window.removeEventListener("popstate",p),window.removeEventListener("beforeunload",r)}return window.addEventListener("popstate",p),window.addEventListener("beforeunload",r,{passive:!0}),{pauseListeners:c,listen:u,destroy:i}}function vp(s,n,a,e=!1,t=!1){return{back:s,current:n,forward:a,replaced:e,position:window.history.length,scroll:t?At():null}}function Lm(s){const{history:n,location:a}=window,e={value:Yr(s,a)},t={value:n.state};t.value||l(e.value,{back:null,current:e.value,forward:null,position:n.length-1,replaced:!0,scroll:null},!0);function l(c,u,r){const i=s.indexOf("#"),d=i>-1?(a.host&&document.querySelector("base")?s:s.slice(i))+c:Am()+s+c;try{n[r?"replaceState":"pushState"](u,"",d),t.value=u}catch(h){console.error(h),a[r?"replace":"assign"](d)}}function o(c,u){const r=Ps({},n.state,vp(t.value.back,c,t.value.forward,!0),u,{position:t.value.position});l(c,r,!0),e.value=c}function p(c,u){const r=Ps({},t.value,n.state,{forward:c,scroll:At()});l(r.current,r,!0);const i=Ps({},vp(e.value,c,null),{position:r.position+1},u);l(c,i,!1),e.value=c}return{location:e,state:t,push:p,replace:o}}function Mm(s){s=Cm(s);const n=Lm(s),a=Rm(s,n.state,n.location,n.replace);function e(l,o=!0){o||a.pauseListeners(),history.go(l)}const t=Ps({location:"",base:s,go:e,createHref:xm.bind(null,s)},n,a);return Object.defineProperty(t,"location",{enumerable:!0,get:()=>n.location.value}),Object.defineProperty(t,"state",{enumerable:!0,get:()=>n.state.value}),t}function Dm(s){return typeof s=="string"||s&&typeof s=="object"}function Qr(s){return typeof s=="string"||typeof s=="symbol"}const Jr=Symbol("");var wp;(function(s){s[s.aborted=4]="aborted",s[s.cancelled=8]="cancelled",s[s.duplicated=16]="duplicated"})(wp||(wp={}));function Wa(s,n){return Ps(new Error,{type:s,[Jr]:!0},n)}function Hn(s,n){return s instanceof Error&&Jr in s&&(n==null||!!(s.type&n))}const Cp="[^/]+?",Om={sensitive:!1,strict:!1,start:!0,end:!0},Im=/[.+*?^${}()[\]/\\]/g;function Nm(s,n){const a=Ps({},Om,n),e=[];let t=a.start?"^":"";const l=[];for(const u of s){const r=u.length?[]:[90];a.strict&&!u.length&&(t+="/");for(let i=0;i<u.length;i++){const d=u[i];let h=40+(a.sensitive?.25:0);if(d.type===0)i||(t+="/"),t+=d.value.replace(Im,"\\$&"),h+=40;else if(d.type===1){const{value:y,repeatable:j,optional:T,regexp:g}=d;l.push({name:y,repeatable:j,optional:T});const m=g||Cp;if(m!==Cp){h+=10;try{new RegExp(`(${m})`)}catch(f){throw new Error(`Invalid custom RegExp for param "${y}" (${m}): `+f.message)}}let w=j?`((?:${m})(?:/(?:${m}))*)`:`(${m})`;i||(w=T&&u.length<2?`(?:/${w})`:"/"+w),T&&(w+="?"),t+=w,h+=20,T&&(h+=-8),j&&(h+=-20),m===".*"&&(h+=-50)}r.push(h)}e.push(r)}if(a.strict&&a.end){const u=e.length-1;e[u][e[u].length-1]+=.7000000000000001}a.strict||(t+="/?"),a.end?t+="$":a.strict&&!t.endsWith("/")&&(t+="(?:/|$)");const o=new RegExp(t,a.sensitive?"":"i");function p(u){const r=u.match(o),i={};if(!r)return null;for(let d=1;d<r.length;d++){const h=r[d]||"",y=l[d-1];i[y.name]=h&&y.repeatable?h.split("/"):h}return i}function c(u){let r="",i=!1;for(const d of s){(!i||!r.endsWith("/"))&&(r+="/"),i=!1;for(const h of d)if(h.type===0)r+=h.value;else if(h.type===1){const{value:y,repeatable:j,optional:T}=h,g=y in u?u[y]:"";if(Mn(g)&&!j)throw new Error(`Provided param "${y}" is an array but it is not repeatable (* or + modifiers)`);const m=Mn(g)?g.join("/"):g;if(!m)if(T)d.length<2&&(r.endsWith("/")?r=r.slice(0,-1):i=!0);else throw new Error(`Missing required param "${y}"`);r+=m}}return r||"/"}return{re:o,score:e,keys:l,parse:p,stringify:c}}function Fm(s,n){let a=0;for(;a<s.length&&a<n.length;){const e=n[a]-s[a];if(e)return e;a++}return s.length<n.length?s.length===1&&s[0]===80?-1:1:s.length>n.length?n.length===1&&n[0]===80?1:-1:0}function Zr(s,n){let a=0;const e=s.score,t=n.score;for(;a<e.length&&a<t.length;){const l=Fm(e[a],t[a]);if(l)return l;a++}if(Math.abs(t.length-e.length)===1){if(kp(e))return 1;if(kp(t))return-1}return t.length-e.length}function kp(s){const n=s[s.length-1];return s.length>0&&n[n.length-1]<0}const $m={type:0,value:""},Bm=/[a-zA-Z0-9_]/;function qm(s){if(!s)return[[]];if(s==="/")return[[$m]];if(!s.startsWith("/"))throw new Error(`Invalid path "${s}"`);function n(h){throw new Error(`ERR (${a})/"${u}": ${h}`)}let a=0,e=a;const t=[];let l;function o(){l&&t.push(l),l=[]}let p=0,c,u="",r="";function i(){u&&(a===0?l.push({type:0,value:u}):a===1||a===2||a===3?(l.length>1&&(c==="*"||c==="+")&&n(`A repeatable param (${u}) must be alone in its segment. eg: '/:ids+.`),l.push({type:1,value:u,regexp:r,repeatable:c==="*"||c==="+",optional:c==="*"||c==="?"})):n("Invalid state to consume buffer"),u="")}function d(){u+=c}for(;p<s.length;){if(c=s[p++],c==="\\"&&a!==2){e=a,a=4;continue}switch(a){case 0:c==="/"?(u&&i(),o()):c===":"?(i(),a=1):d();break;case 4:d(),a=e;break;case 1:c==="("?a=2:Bm.test(c)?d():(i(),a=0,c!=="*"&&c!=="?"&&c!=="+"&&p--);break;case 2:c===")"?r[r.length-1]=="\\"?r=r.slice(0,-1)+c:a=3:r+=c;break;case 3:i(),a=0,c!=="*"&&c!=="?"&&c!=="+"&&p--,r="";break;default:n("Unknown state");break}}return a===2&&n(`Unfinished custom RegExp for param "${u}"`),i(),o(),t}function zm(s,n,a){const e=Nm(qm(s.path),a),t=Ps(e,{record:s,parent:n,children:[],alias:[]});return n&&!t.record.aliasOf==!n.record.aliasOf&&n.children.push(t),t}function Vm(s,n){const a=[],e=new Map;n=Sp({strict:!1,end:!0,sensitive:!1},n);function t(i){return e.get(i)}function l(i,d,h){const y=!h,j=Pp(i);j.aliasOf=h&&h.record;const T=Sp(n,i),g=[j];if("alias"in i){const f=typeof i.alias=="string"?[i.alias]:i.alias;for(const x of f)g.push(Pp(Ps({},j,{components:h?h.record.components:j.components,path:x,aliasOf:h?h.record:j})))}let m,w;for(const f of g){const{path:x}=f;if(d&&x[0]!=="/"){const k=d.record.path,_=k[k.length-1]==="/"?"":"/";f.path=d.record.path+(x&&_+x)}if(m=zm(f,d,T),h?h.alias.push(m):(w=w||m,w!==m&&w.alias.push(m),y&&i.name&&!Ep(m)&&o(i.name)),si(m)&&c(m),j.children){const k=j.children;for(let _=0;_<k.length;_++)l(k[_],m,h&&h.children[_])}h=h||m}return w?()=>{o(w)}:he}function o(i){if(Qr(i)){const d=e.get(i);d&&(e.delete(i),a.splice(a.indexOf(d),1),d.children.forEach(o),d.alias.forEach(o))}else{const d=a.indexOf(i);d>-1&&(a.splice(d,1),i.record.name&&e.delete(i.record.name),i.children.forEach(o),i.alias.forEach(o))}}function p(){return a}function c(i){const d=Gm(i,a);a.splice(d,0,i),i.record.name&&!Ep(i)&&e.set(i.record.name,i)}function u(i,d){let h,y={},j,T;if("name"in i&&i.name){if(h=e.get(i.name),!h)throw Wa(1,{location:i});T=h.record.name,y=Ps(xp(d.params,h.keys.filter(w=>!w.optional).concat(h.parent?h.parent.keys.filter(w=>w.optional):[]).map(w=>w.name)),i.params&&xp(i.params,h.keys.map(w=>w.name))),j=h.stringify(y)}else if(i.path!=null)j=i.path,h=a.find(w=>w.re.test(j)),h&&(y=h.parse(j),T=h.record.name);else{if(h=d.name?e.get(d.name):a.find(w=>w.re.test(d.path)),!h)throw Wa(1,{location:i,currentLocation:d});T=h.record.name,y=Ps({},d.params,i.params),j=h.stringify(y)}const g=[];let m=h;for(;m;)g.unshift(m.record),m=m.parent;return{name:T,path:j,params:y,matched:g,meta:Hm(g)}}s.forEach(i=>l(i));function r(){a.length=0,e.clear()}return{addRoute:l,resolve:u,removeRoute:o,clearRoutes:r,getRoutes:p,getRecordMatcher:t}}function xp(s,n){const a={};for(const e of n)e in s&&(a[e]=s[e]);return a}function Pp(s){const n={path:s.path,redirect:s.redirect,name:s.name,meta:s.meta||{},aliasOf:s.aliasOf,beforeEnter:s.beforeEnter,props:Um(s),children:s.children||[],instances:{},leaveGuards:new Set,updateGuards:new Set,enterCallbacks:{},components:"components"in s?s.components||null:s.component&&{default:s.component}};return Object.defineProperty(n,"mods",{value:{}}),n}function Um(s){const n={},a=s.props||!1;if("component"in s)n.default=a;else for(const e in s.components)n[e]=typeof a=="object"?a[e]:a;return n}function Ep(s){for(;s;){if(s.record.aliasOf)return!0;s=s.parent}return!1}function Hm(s){return s.reduce((n,a)=>Ps(n,a.meta),{})}function Sp(s,n){const a={};for(const e in s)a[e]=e in n?n[e]:s[e];return a}function Gm(s,n){let a=0,e=n.length;for(;a!==e;){const l=a+e>>1;Zr(s,n[l])<0?e=l:a=l+1}const t=Wm(s);return t&&(e=n.lastIndexOf(t,e-1)),e}function Wm(s){let n=s;for(;n=n.parent;)if(si(n)&&Zr(s,n)===0)return n}function si({record:s}){return!!(s.name||s.components&&Object.keys(s.components).length||s.redirect)}function Km(s){const n={};if(s===""||s==="?")return n;const e=(s[0]==="?"?s.slice(1):s).split("&");for(let t=0;t<e.length;++t){const l=e[t].replace(Hr," "),o=l.indexOf("="),p=we(o<0?l:l.slice(0,o)),c=o<0?null:we(l.slice(o+1));if(p in n){let u=n[p];Mn(u)||(u=n[p]=[u]),u.push(c)}else n[p]=c}return n}function Tp(s){let n="";for(let a in s){const e=s[a];if(a=mm(a),e==null){e!==void 0&&(n+=(n.length?"&":"")+a);continue}(Mn(e)?e.map(l=>l&&wl(l)):[e&&wl(e)]).forEach(l=>{l!==void 0&&(n+=(n.length?"&":"")+a,l!=null&&(n+="="+l))})}return n}function Xm(s){const n={};for(const a in s){const e=s[a];e!==void 0&&(n[a]=Mn(e)?e.map(t=>t==null?null:""+t):e==null?e:""+e)}return n}const Ym=Symbol(""),Ap=Symbol(""),Rt=Symbol(""),mo=Symbol(""),kl=Symbol("");function se(){let s=[];function n(e){return s.push(e),()=>{const t=s.indexOf(e);t>-1&&s.splice(t,1)}}function a(){s=[]}return{add:n,list:()=>s.slice(),reset:a}}function ia(s,n,a,e,t,l=o=>o()){const o=e&&(e.enterCallbacks[t]=e.enterCallbacks[t]||[]);return()=>new Promise((p,c)=>{const u=d=>{d===!1?c(Wa(4,{from:a,to:n})):d instanceof Error?c(d):Dm(d)?c(Wa(2,{from:n,to:d})):(o&&e.enterCallbacks[t]===o&&typeof d=="function"&&o.push(d),p())},r=l(()=>s.call(e&&e.instances[t],n,a,u));let i=Promise.resolve(r);s.length<3&&(i=i.then(u)),i.catch(d=>c(d))})}function sl(s,n,a,e,t=l=>l()){const l=[];for(const o of s)for(const p in o.components){let c=o.components[p];if(!(n!=="beforeRouteEnter"&&!o.instances[p]))if(Vr(c)){const r=(c.__vccOpts||c)[n];r&&l.push(ia(r,a,e,o,p,t))}else{let u=c();l.push(()=>u.then(r=>{if(!r)throw new Error(`Couldn't resolve component "${p}" at "${o.path}"`);const i=em(r)?r.default:r;o.mods[p]=r,o.components[p]=i;const h=(i.__vccOpts||i)[n];return h&&ia(h,a,e,o,p,t)()}))}}return l}function Rp(s){const n=jn(Rt),a=jn(mo),e=ss(()=>{const c=cs(s.to);return n.resolve(c)}),t=ss(()=>{const{matched:c}=e.value,{length:u}=c,r=c[u-1],i=a.matched;if(!r||!i.length)return-1;const d=i.findIndex(Ga.bind(null,r));if(d>-1)return d;const h=Lp(c[u-2]);return u>1&&Lp(r)===h&&i[i.length-1].path!==h?i.findIndex(Ga.bind(null,c[u-2])):d}),l=ss(()=>t.value>-1&&nj(a.params,e.value.params)),o=ss(()=>t.value>-1&&t.value===a.matched.length-1&&Xr(a.params,e.value.params));function p(c={}){if(sj(c)){const u=n[cs(s.replace)?"replace":"push"](cs(s.to)).catch(he);return s.viewTransition&&typeof document<"u"&&"startViewTransition"in document&&document.startViewTransition(()=>u),u}return Promise.resolve()}return{route:e,href:ss(()=>e.value.href),isActive:l,isExactActive:o,navigate:p}}function Qm(s){return s.length===1?s[0]:s}const Jm=Gs({name:"RouterLink",compatConfig:{MODE:3},props:{to:{type:[String,Object],required:!0},replace:Boolean,activeClass:String,exactActiveClass:String,custom:Boolean,ariaCurrentValue:{type:String,default:"page"},viewTransition:Boolean},useLink:Rp,setup(s,{slots:n}){const a=Se(Rp(s)),{options:e}=jn(Rt),t=ss(()=>({[Mp(s.activeClass,e.linkActiveClass,"router-link-active")]:a.isActive,[Mp(s.exactActiveClass,e.linkExactActiveClass,"router-link-exact-active")]:a.isExactActive}));return()=>{const l=n.default&&Qm(n.default(a));return s.custom?l:Oe("a",{"aria-current":a.isExactActive?s.ariaCurrentValue:null,href:a.href,onClick:a.navigate,class:t.value},l)}}}),Zm=Jm;function sj(s){if(!(s.metaKey||s.altKey||s.ctrlKey||s.shiftKey)&&!s.defaultPrevented&&!(s.button!==void 0&&s.button!==0)){if(s.currentTarget&&s.currentTarget.getAttribute){const n=s.currentTarget.getAttribute("target");if(/\b_blank\b/i.test(n))return}return s.preventDefault&&s.preventDefault(),!0}}function nj(s,n){for(const a in n){const e=n[a],t=s[a];if(typeof e=="string"){if(e!==t)return!1}else if(!Mn(t)||t.length!==e.length||e.some((l,o)=>l!==t[o]))return!1}return!0}function Lp(s){return s?s.aliasOf?s.aliasOf.path:s.path:""}const Mp=(s,n,a)=>s??n??a,aj=Gs({name:"RouterView",inheritAttrs:!1,props:{name:{type:String,default:"default"},route:Object},compatConfig:{MODE:3},setup(s,{attrs:n,slots:a}){const e=jn(kl),t=ss(()=>s.route||e.value),l=jn(Ap,0),o=ss(()=>{let u=cs(l);const{matched:r}=t.value;let i;for(;(i=r[u])&&!i.components;)u++;return u}),p=ss(()=>t.value.matched[o.value]);Ze(Ap,ss(()=>o.value+1)),Ze(Ym,p),Ze(kl,t);const c=os();return Hs(()=>[c.value,p.value,s.name],([u,r,i],[d,h,y])=>{r&&(r.instances[i]=u,h&&h!==r&&u&&u===d&&(r.leaveGuards.size||(r.leaveGuards=h.leaveGuards),r.updateGuards.size||(r.updateGuards=h.updateGuards))),u&&r&&(!h||!Ga(r,h)||!d)&&(r.enterCallbacks[i]||[]).forEach(j=>j(u))},{flush:"post"}),()=>{const u=t.value,r=s.name,i=p.value,d=i&&i.components[r];if(!d)return Dp(a.default,{Component:d,route:u});const h=i.props[r],y=h?h===!0?u.params:typeof h=="function"?h(u):h:null,T=Oe(d,Ps({},y,n,{onVnodeUnmounted:g=>{g.component.isUnmounted&&(i.instances[r]=null)},ref:c}));return Dp(a.default,{Component:T,route:u})||T}}});function Dp(s,n){if(!s)return null;const a=s(n);return a.length===1?a[0]:a}const ej=aj;function tj(s){const n=Vm(s.routes,s),a=s.parseQuery||Km,e=s.stringifyQuery||Tp,t=s.history,l=se(),o=se(),p=se(),c=so(la);let u=la;Fa&&s.scrollBehavior&&"scrollRestoration"in history&&(history.scrollRestoration="manual");const r=Jt.bind(null,O=>""+O),i=Jt.bind(null,fm),d=Jt.bind(null,we);function h(O,W){let H,Q;return Qr(O)?(H=n.getRecordMatcher(O),Q=W):Q=O,n.addRoute(Q,H)}function y(O){const W=n.getRecordMatcher(O);W&&n.removeRoute(W)}function j(){return n.getRoutes().map(O=>O.record)}function T(O){return!!n.getRecordMatcher(O)}function g(O,W){if(W=Ps({},W||c.value),typeof O=="string"){const R=Zt(a,O,W.path),$=n.resolve({path:R.path},W),V=t.createHref(R.fullPath);return Ps(R,$,{params:d($.params),hash:we(R.hash),redirectedFrom:void 0,href:V})}let H;if(O.path!=null)H=Ps({},O,{path:Zt(a,O.path,W.path).path});else{const R=Ps({},O.params);for(const $ in R)R[$]==null&&delete R[$];H=Ps({},O,{params:i(R)}),W.params=i(W.params)}const Q=n.resolve(H,W),gs=O.hash||"";Q.params=r(d(Q.params));const C=_m(e,Ps({},O,{hash:dm(gs),path:Q.path})),P=t.createHref(C);return Ps({fullPath:C,hash:gs,query:e===Tp?Xm(O.query):O.query||{}},Q,{redirectedFrom:void 0,href:P})}function m(O){return typeof O=="string"?Zt(a,O,c.value.path):Ps({},O)}function w(O,W){if(u!==O)return Wa(8,{from:W,to:O})}function f(O){return _(O)}function x(O){return f(Ps(m(O),{replace:!0}))}function k(O){const W=O.matched[O.matched.length-1];if(W&&W.redirect){const{redirect:H}=W;let Q=typeof H=="function"?H(O):H;return typeof Q=="string"&&(Q=Q.includes("?")||Q.includes("#")?Q=m(Q):{path:Q},Q.params={}),Ps({query:O.query,hash:O.hash,params:Q.path!=null?{}:O.params},Q)}}function _(O,W){const H=u=g(O),Q=c.value,gs=O.state,C=O.force,P=O.replace===!0,R=k(H);if(R)return _(Ps(m(R),{state:typeof R=="object"?Ps({},gs,R.state):gs,force:C,replace:P}),W||H);const $=H;$.redirectedFrom=W;let V;return!C&&ym(e,Q,H)&&(V=Wa(16,{to:$,from:Q}),Rs(Q,Q,!0,!1)),(V?Promise.resolve(V):I($,Q)).catch(q=>Hn(q)?Hn(q,2)?q:Xs(q):us(q,$,Q)).then(q=>{if(q){if(Hn(q,2))return _(Ps({replace:P},m(q.to),{state:typeof q.to=="object"?Ps({},gs,q.to.state):gs,force:C}),W||$)}else q=U($,Q,!0,P,gs);return J($,Q,q),q})}function D(O,W){const H=w(O,W);return H?Promise.reject(H):Promise.resolve()}function S(O){const W=rs.values().next().value;return W&&typeof W.runWithContext=="function"?W.runWithContext(O):O()}function I(O,W){let H;const[Q,gs,C]=lj(O,W);H=sl(Q.reverse(),"beforeRouteLeave",O,W);for(const R of Q)R.leaveGuards.forEach($=>{H.push(ia($,O,W))});const P=D.bind(null,O,W);return H.push(P),_s(H).then(()=>{H=[];for(const R of l.list())H.push(ia(R,O,W));return H.push(P),_s(H)}).then(()=>{H=sl(gs,"beforeRouteUpdate",O,W);for(const R of gs)R.updateGuards.forEach($=>{H.push(ia($,O,W))});return H.push(P),_s(H)}).then(()=>{H=[];for(const R of C)if(R.beforeEnter)if(Mn(R.beforeEnter))for(const $ of R.beforeEnter)H.push(ia($,O,W));else H.push(ia(R.beforeEnter,O,W));return H.push(P),_s(H)}).then(()=>(O.matched.forEach(R=>R.enterCallbacks={}),H=sl(C,"beforeRouteEnter",O,W,S),H.push(P),_s(H))).then(()=>{H=[];for(const R of o.list())H.push(ia(R,O,W));return H.push(P),_s(H)}).catch(R=>Hn(R,8)?R:Promise.reject(R))}function J(O,W,H){p.list().forEach(Q=>S(()=>Q(O,W,H)))}function U(O,W,H,Q,gs){const C=w(O,W);if(C)return C;const P=W===la,R=Fa?history.state:{};H&&(Q||P?t.replace(O.fullPath,Ps({scroll:P&&R&&R.scroll},gs)):t.push(O.fullPath,gs)),c.value=O,Rs(O,W,H,P),Xs()}let F;function G(){F||(F=t.listen((O,W,H)=>{if(!fs.listening)return;const Q=g(O),gs=k(Q);if(gs){_(Ps(gs,{replace:!0,force:!0}),Q).catch(he);return}u=Q;const C=c.value;Fa&&Sm(yp(C.fullPath,H.delta),At()),I(Q,C).catch(P=>Hn(P,12)?P:Hn(P,2)?(_(Ps(m(P.to),{force:!0}),Q).then(R=>{Hn(R,20)&&!H.delta&&H.type===Ce.pop&&t.go(-1,!1)}).catch(he),Promise.reject()):(H.delta&&t.go(-H.delta,!1),us(P,Q,C))).then(P=>{P=P||U(Q,C,!1),P&&(H.delta&&!Hn(P,8)?t.go(-H.delta,!1):H.type===Ce.pop&&Hn(P,20)&&t.go(-1,!1)),J(Q,C,P)}).catch(he)}))}let is=se(),ps=se(),Y;function us(O,W,H){Xs(O);const Q=ps.list();return Q.length?Q.forEach(gs=>gs(O,W,H)):console.error(O),Promise.reject(O)}function Ns(){return Y&&c.value!==la?Promise.resolve():new Promise((O,W)=>{is.add([O,W])})}function Xs(O){return Y||(Y=!O,G(),is.list().forEach(([W,H])=>O?H(O):W()),is.reset()),O}function Rs(O,W,H,Q){const{scrollBehavior:gs}=s;if(!Fa||!gs)return Promise.resolve();const C=!H&&Tm(yp(O.fullPath,0))||(Q||!H)&&history.state&&history.state.scroll||null;return Js().then(()=>gs(O,W,C)).then(P=>P&&Em(P)).catch(P=>us(P,O,W))}const zs=O=>t.go(O);let ts;const rs=new Set,fs={currentRoute:c,listening:!0,addRoute:h,removeRoute:y,clearRoutes:n.clearRoutes,hasRoute:T,getRoutes:j,resolve:g,options:s,push:f,replace:x,go:zs,back:()=>zs(-1),forward:()=>zs(1),beforeEach:l.add,beforeResolve:o.add,afterEach:p.add,onError:ps.add,isReady:Ns,install(O){const W=this;O.component("RouterLink",Zm),O.component("RouterView",ej),O.config.globalProperties.$router=W,Object.defineProperty(O.config.globalProperties,"$route",{enumerable:!0,get:()=>cs(c)}),Fa&&!ts&&c.value===la&&(ts=!0,f(t.location).catch(gs=>{}));const H={};for(const gs in la)Object.defineProperty(H,gs,{get:()=>c.value[gs],enumerable:!0});O.provide(Rt,W),O.provide(mo,Qc(H)),O.provide(kl,c);const Q=O.unmount;rs.add(O),O.unmount=function(){rs.delete(O),rs.size<1&&(u=la,F&&F(),F=null,c.value=la,ts=!1,Y=!1),Q()}}};function _s(O){return O.reduce((W,H)=>W.then(()=>S(H)),Promise.resolve())}return fs}function lj(s,n){const a=[],e=[],t=[],l=Math.max(n.matched.length,s.matched.length);for(let o=0;o<l;o++){const p=n.matched[o];p&&(s.matched.find(u=>Ga(u,p))?e.push(p):a.push(p));const c=s.matched[o];c&&(n.matched.find(u=>Ga(u,c))||t.push(c))}return[a,e,t]}function Ie(){return jn(Rt)}function Ka(s){return jn(mo)}function oj(s){return document.readyState==="loading"?new Promise(n=>{document.addEventListener("DOMContentLoaded",()=>n(s))}):Promise.resolve(s)}const pj=Gs({setup(s,{slots:n}){const a=os(!1);return wn(()=>a.value=!0),()=>a.value?n.default&&n.default({}):n.placeholder&&n.placeholder({})}});function cj(s){try{return JSON.parse(s||"{}")}catch(n){return console.error("[SSG] On state deserialization -",n,s),{}}}function rj(s,n,a,e){const{transformState:t,registerComponents:l=!0,useHead:o=!0,rootContainer:p="#app"}={};async function c(u){const r=Kd(s);let i;o&&r.use(i=am());const d=tj({history:Mm(n.base),...n}),{routes:h}=n;l&&r.component("ClientOnly",pj);const y=[],g={app:r,head:i,isClient:!0,router:d,routes:h,onSSRAppRendered:()=>{},triggerOnSSRAppRendered:()=>Promise.all(y.map(x=>x())),initialState:{},transformState:t,routePath:u};await oj(),g.initialState=t?.(window.__INITIAL_STATE__||{})||cj(window.__INITIAL_STATE__),await a?.(g),r.use(d);let m,w=!0;d.beforeEach((x,k,_)=>{(w||m&&m===x.path)&&(w=!1,m=x.path,x.meta.state=g.initialState),_()});const f=g.initialState;return{...g,initialState:f}}return(async()=>{const{app:u,router:r}=await c();await r.isReady(),u.mount(p,!0)})(),c}/*!
 * pinia v3.0.3
 * (c) 2025 Eduardo San Martin Morote
 * @license MIT
 */let ni;const Lt=s=>ni=s,ai=Symbol();function xl(s){return s&&typeof s=="object"&&Object.prototype.toString.call(s)==="[object Object]"&&typeof s.toJSON!="function"}var me;(function(s){s.direct="direct",s.patchObject="patch object",s.patchFunction="patch function"})(me||(me={}));function ij(){const s=Hl(!0),n=s.run(()=>os({}));let a=[],e=[];const t=Zl({install(l){Lt(t),t._a=l,l.provide(ai,t),l.config.globalProperties.$pinia=t,e.forEach(o=>a.push(o)),e=[]},use(l){return this._a?a.push(l):e.push(l),this},_p:a,_a:null,_e:s,_s:new Map,state:n});return t}const ei=()=>{};function Op(s,n,a,e=ei){s.push(n);const t=()=>{const l=s.indexOf(n);l>-1&&(s.splice(l,1),e())};return!a&&Ic()&&Mu(t),t}function Da(s,...n){s.slice().forEach(a=>{a(...n)})}const uj=s=>s(),Ip=Symbol(),nl=Symbol();function Pl(s,n){s instanceof Map&&n instanceof Map?n.forEach((a,e)=>s.set(e,a)):s instanceof Set&&n instanceof Set&&n.forEach(s.add,s);for(const a in n){if(!n.hasOwnProperty(a))continue;const e=n[a],t=s[a];xl(t)&&xl(e)&&s.hasOwnProperty(a)&&!Fs(e)&&!da(e)?s[a]=Pl(t,e):s[a]=e}return s}const hj=Symbol();function dj(s){return!xl(s)||!Object.prototype.hasOwnProperty.call(s,hj)}const{assign:pa}=Object;function mj(s){return!!(Fs(s)&&s.effect)}function jj(s,n,a,e){const{state:t,actions:l,getters:o}=n,p=a.state.value[s];let c;function u(){p||(a.state.value[s]=t?t():{});const r=ah(a.state.value[s]);return pa(r,l,Object.keys(o||{}).reduce((i,d)=>(i[d]=Zl(ss(()=>{Lt(a);const h=a._s.get(s);return o[d].call(h,h)})),i),{}))}return c=ti(s,u,n,a,e,!0),c}function ti(s,n,a={},e,t,l){let o;const p=pa({actions:{}},a),c={deep:!0};let u,r,i=[],d=[],h;const y=e.state.value[s];!l&&!y&&(e.state.value[s]={}),os({});let j;function T(D){let S;u=r=!1,typeof D=="function"?(D(e.state.value[s]),S={type:me.patchFunction,storeId:s,events:h}):(Pl(e.state.value[s],D),S={type:me.patchObject,payload:D,storeId:s,events:h});const I=j=Symbol();Js().then(()=>{j===I&&(u=!0)}),r=!0,Da(i,S,e.state.value[s])}const g=l?function(){const{state:S}=a,I=S?S():{};this.$patch(J=>{pa(J,I)})}:ei;function m(){o.stop(),i=[],d=[],e._s.delete(s)}const w=(D,S="")=>{if(Ip in D)return D[nl]=S,D;const I=function(){Lt(e);const J=Array.from(arguments),U=[],F=[];function G(Y){U.push(Y)}function is(Y){F.push(Y)}Da(d,{args:J,name:I[nl],store:x,after:G,onError:is});let ps;try{ps=D.apply(this&&this.$id===s?this:x,J)}catch(Y){throw Da(F,Y),Y}return ps instanceof Promise?ps.then(Y=>(Da(U,Y),Y)).catch(Y=>(Da(F,Y),Promise.reject(Y))):(Da(U,ps),ps)};return I[Ip]=!0,I[nl]=S,I},f={_p:e,$id:s,$onAction:Op.bind(null,d),$patch:T,$reset:g,$subscribe(D,S={}){const I=Op(i,D,S.detached,()=>J()),J=o.run(()=>Hs(()=>e.state.value[s],U=>{(S.flush==="sync"?r:u)&&D({storeId:s,type:me.direct,events:h},U)},pa({},c,S)));return I},$dispose:m},x=Se(f);e._s.set(s,x);const _=(e._a&&e._a.runWithContext||uj)(()=>e._e.run(()=>(o=Hl()).run(()=>n({action:w}))));for(const D in _){const S=_[D];if(Fs(S)&&!mj(S)||da(S))l||(y&&dj(S)&&(Fs(S)?S.value=y[D]:Pl(S,y[D])),e.state.value[s][D]=S);else if(typeof S=="function"){const I=w(S,D);_[D]=I,p.actions[D]=S}}return pa(x,_),pa(xs(x),_),Object.defineProperty(x,"$state",{get:()=>e.state.value[s],set:D=>{T(S=>{pa(S,D)})}}),e._p.forEach(D=>{pa(x,o.run(()=>D({store:x,app:e._a,pinia:e,options:p})))}),y&&l&&a.hydrate&&a.hydrate(x.$state,y),u=!0,r=!0,x}/*! #__NO_SIDE_EFFECTS__ */function fj(s,n,a){let e;const t=typeof n=="function";e=t?a:n;function l(o,p){const c=vr();return o=o||(c?jn(ai,null):null),o&&Lt(o),o=ni,o._s.has(s)||(t?ti(s,n,e,o):jj(s,e,o)),o._s.get(s)}return l.$id=s,l}/*!
  * vue-i18n v12.0.0-alpha.3
  * (c) 2016-present kazuya kawaguchi and contributors
  * Released under the MIT License.
  */const na=typeof window<"u";let yn,Ra;{const s=na&&window.performance;s&&s.mark&&s.measure&&s.clearMarks&&s.clearMeasures&&(yn=n=>{s.mark(n)},Ra=(n,a,e)=>{s.measure(n,a,e),s.clearMarks(a),s.clearMarks(e)})}const gj=/\{([0-9a-zA-Z]+)\}/g;function Mt(s,...n){return n.length===1&&ks(n[0])&&(n=n[0]),(!n||!n.hasOwnProperty)&&(n={}),s.replace(gj,(a,e)=>n.hasOwnProperty(e)?n[e]:"")}const zn=(s,n=!1)=>n?Symbol.for(s):Symbol(s),bj=(s,n,a)=>_j({l:s,k:n,s:a}),_j=s=>JSON.stringify(s).replace(/\u2028/g,"\\u2028").replace(/\u2029/g,"\\u2029").replace(/\u0027/g,"\\u0027"),Ws=s=>typeof s=="number"&&isFinite(s),yj=s=>jo(s)==="[object Date]",ft=s=>jo(s)==="[object RegExp]",Dt=s=>Cs(s)&&Object.keys(s).length===0,sn=Object.assign,vj=Object.create,Ts=(s=null)=>vj(s);let Np;const wj=()=>Np||(Np=typeof globalThis<"u"?globalThis:typeof self<"u"?self:typeof window<"u"?window:typeof global<"u"?global:Ts());function Fp(s){return s.replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;").replace(/'/g,"&apos;")}const Cj=Object.prototype.hasOwnProperty;function An(s,n){return Cj.call(s,n)}const qs=Array.isArray,Os=s=>typeof s=="function",es=s=>typeof s=="string",Is=s=>typeof s=="boolean",ks=s=>s!==null&&typeof s=="object",kj=s=>ks(s)&&Os(s.then)&&Os(s.catch),li=Object.prototype.toString,jo=s=>li.call(s),Cs=s=>jo(s)==="[object Object]",xj=s=>s==null?"":qs(s)||Cs(s)&&s.toString===li?JSON.stringify(s,null,2):String(s);function oi(s,n=""){return s.reduce((a,e,t)=>t===0?a+e:a+n+e,"")}const $p=2;function Pj(s,n=0,a=s.length){const e=s.split(/\r?\n/);let t=0;const l=[];for(let o=0;o<e.length;o++)if(t+=e[o].length+1,t>=n){for(let p=o-$p;p<=o+$p||a>t;p++){if(p<0||p>=e.length)continue;const c=p+1;l.push(`${c}${" ".repeat(3-String(c).length)}|  ${e[p]}`);const u=e[p].length;if(p===o){const r=n-(t-u)+1,i=Math.max(1,a>t?u-r:a-n);l.push("   |  "+" ".repeat(r)+"^".repeat(i))}else if(p>o){if(a>t){const r=Math.max(Math.min(a-t,u),1);l.push("   |  "+"^".repeat(r))}t+=u+1}}break}return l.join(`
`)}function _a(s,n){typeof console<"u"&&(console.warn("[intlify] "+s),n&&console.warn(n.stack))}const Bp={};function Ej(s){Bp[s]||(Bp[s]=!0,_a(s))}function pi(){const s=new Map;return{events:s,on(a,e){const t=s.get(a);t&&t.push(e)||s.set(a,[e])},off(a,e){const t=s.get(a);t&&t.splice(t.indexOf(e)>>>0,1)},emit(a,e){(s.get(a)||[]).slice().map(t=>t(e)),(s.get("*")||[]).slice().map(t=>t(a,e))}}}const Ge=s=>!ks(s)||qs(s);function et(s,n){if(Ge(s)||Ge(n))throw new Error("Invalid value");const a=[{src:s,des:n}];for(;a.length;){const{src:e,des:t}=a.pop();Object.keys(e).forEach(l=>{l!=="__proto__"&&(ks(e[l])&&!ks(t[l])&&(t[l]=Array.isArray(e[l])?[]:Ts()),Ge(t[l])||Ge(e[l])?t[l]=e[l]:a.push({src:e[l],des:t[l]}))})}}function Sj(s,n,a){return{line:s,column:n,offset:a}}function El(s,n,a){return{start:s,end:n}}const ms={EXPECTED_TOKEN:1,INVALID_TOKEN_IN_PLACEHOLDER:2,UNTERMINATED_SINGLE_QUOTE_IN_PLACEHOLDER:3,UNKNOWN_ESCAPE_SEQUENCE:4,INVALID_UNICODE_ESCAPE_SEQUENCE:5,UNBALANCED_CLOSING_BRACE:6,UNTERMINATED_CLOSING_BRACE:7,EMPTY_PLACEHOLDER:8,NOT_ALLOW_NEST_PLACEHOLDER:9,INVALID_LINKED_FORMAT:10,MUST_HAVE_MESSAGES_IN_PLURAL:11,UNEXPECTED_EMPTY_LINKED_MODIFIER:12,UNEXPECTED_EMPTY_LINKED_KEY:13,UNEXPECTED_LEXICAL_ANALYSIS:14,UNHANDLED_CODEGEN_NODE_TYPE:15,UNHANDLED_MINIFIER_NODE_TYPE:16},Tj=17,Aj={[ms.EXPECTED_TOKEN]:"Expected token: '{0}'",[ms.INVALID_TOKEN_IN_PLACEHOLDER]:"Invalid token in placeholder: '{0}'",[ms.UNTERMINATED_SINGLE_QUOTE_IN_PLACEHOLDER]:"Unterminated single quote in placeholder",[ms.UNKNOWN_ESCAPE_SEQUENCE]:"Unknown escape sequence: \\{0}",[ms.INVALID_UNICODE_ESCAPE_SEQUENCE]:"Invalid unicode escape sequence: {0}",[ms.UNBALANCED_CLOSING_BRACE]:"Unbalanced closing brace",[ms.UNTERMINATED_CLOSING_BRACE]:"Unterminated closing brace",[ms.EMPTY_PLACEHOLDER]:"Empty placeholder",[ms.NOT_ALLOW_NEST_PLACEHOLDER]:"Not allowed nest placeholder",[ms.INVALID_LINKED_FORMAT]:"Invalid linked format",[ms.MUST_HAVE_MESSAGES_IN_PLURAL]:"Plural must have messages",[ms.UNEXPECTED_EMPTY_LINKED_MODIFIER]:"Unexpected empty linked modifier",[ms.UNEXPECTED_EMPTY_LINKED_KEY]:"Unexpected empty linked key",[ms.UNEXPECTED_LEXICAL_ANALYSIS]:"Unexpected lexical analysis in token: '{0}'",[ms.UNHANDLED_CODEGEN_NODE_TYPE]:"unhandled codegen node type: '{0}'",[ms.UNHANDLED_MINIFIER_NODE_TYPE]:"unhandled mimifier node type: '{0}'"};function Ne(s,n,a={}){const{domain:e,messages:t,args:l}=a,o=Mt((t||Aj)[s]||"",...l||[]),p=new SyntaxError(String(o));return p.code=s,n&&(p.location=n),p.domain=e,p}function Rj(s){throw s}const Lj="minifier";function $a(s){switch(s.t=s.type,s.type){case 0:{const n=s;$a(n.body),n.b=n.body,delete n.body;break}case 1:{const n=s,a=n.cases;for(let e=0;e<a.length;e++)$a(a[e]);n.c=a,delete n.cases;break}case 2:{const n=s,a=n.items;for(let e=0;e<a.length;e++)$a(a[e]);n.i=a,delete n.items,n.static&&(n.s=n.static,delete n.static);break}case 3:case 9:case 8:case 7:{const n=s;n.value&&(n.v=n.value,delete n.value);break}case 6:{const n=s;$a(n.key),n.k=n.key,delete n.key,n.modifier&&($a(n.modifier),n.m=n.modifier,delete n.modifier);break}case 5:{const n=s;n.i=n.index,delete n.index;break}case 4:{const n=s;n.k=n.key,delete n.key;break}default:throw Ne(ms.UNHANDLED_MINIFIER_NODE_TYPE,null,{domain:Lj,args:[s.type]})}delete s.type}function Mj(s){const n=s.body;return n.type===2?qp(n):n.cases.forEach(a=>qp(a)),s}function qp(s){if(s.items.length===1){const n=s.items[0];(n.type===3||n.type===9)&&(s.static=n.value,delete n.value)}else{const n=[];for(let a=0;a<s.items.length;a++){const e=s.items[a];if(!(e.type===3||e.type===9)||e.value==null)break;n.push(e.value)}if(n.length===s.items.length){s.static=oi(n);for(let a=0;a<s.items.length;a++){const e=s.items[a];(e.type===3||e.type===9)&&delete e.value}}}}const Gn=" ",Dj="\r",cn=`
`,Oj="\u2028",Ij="\u2029";function Nj(s){const n=s;let a=0,e=1,t=1,l=0;const o=_=>n[_]===Dj&&n[_+1]===cn,p=_=>n[_]===cn,c=_=>n[_]===Ij,u=_=>n[_]===Oj,r=_=>o(_)||p(_)||c(_)||u(_),i=()=>a,d=()=>e,h=()=>t,y=()=>l,j=_=>o(_)||c(_)||u(_)?cn:n[_],T=()=>j(a),g=()=>j(a+l);function m(){return l=0,r(a)&&(e++,t=0),o(a)&&a++,a++,t++,n[a]}function w(){return o(a+l)&&l++,l++,n[a+l]}function f(){a=0,e=1,t=1,l=0}function x(_=0){l=_}function k(){const _=a+l;for(;_!==a;)m();l=0}return{index:i,line:d,column:h,peekOffset:y,charAt:j,currentChar:T,currentPeek:g,next:m,peek:w,reset:f,resetPeek:x,skipToPeek:k}}const oa=void 0,Fj=".",zp="'",$j="tokenizer";function Bj(s,n={}){const a=n.location!==!1,e=Nj(s),t=()=>e.index(),l=()=>Sj(e.line(),e.column(),e.index()),o=l(),p=t(),c={currentType:13,offset:p,startLoc:o,endLoc:o,lastType:13,lastOffset:p,lastStartLoc:o,lastEndLoc:o,braceNest:0,inLinked:!1,text:""},u=()=>c,{onError:r}=n;function i(b,v,A,...M){const Z=u();if(v.column+=A,v.offset+=A,r){const K=a?El(Z.startLoc,v):null,ns=Ne(b,K,{domain:$j,args:M});r(ns)}}function d(b,v,A){b.endLoc=l(),b.currentType=v;const M={type:v};return a&&(M.loc=El(b.startLoc,b.endLoc)),A!=null&&(M.value=A),M}const h=b=>d(b,13);function y(b,v){return b.currentChar()===v?(b.next(),v):(i(ms.EXPECTED_TOKEN,l(),0,v),"")}function j(b){let v="";for(;b.currentPeek()===Gn||b.currentPeek()===cn;)v+=b.currentPeek(),b.peek();return v}function T(b){const v=j(b);return b.skipToPeek(),v}function g(b){if(b===oa)return!1;const v=b.charCodeAt(0);return v>=97&&v<=122||v>=65&&v<=90||v===95}function m(b){if(b===oa)return!1;const v=b.charCodeAt(0);return v>=48&&v<=57}function w(b,v){const{currentType:A}=v;if(A!==2)return!1;j(b);const M=g(b.currentPeek());return b.resetPeek(),M}function f(b,v){const{currentType:A}=v;if(A!==2)return!1;j(b);const M=b.currentPeek()==="-"?b.peek():b.currentPeek(),Z=m(M);return b.resetPeek(),Z}function x(b,v){const{currentType:A}=v;if(A!==2)return!1;j(b);const M=b.currentPeek()===zp;return b.resetPeek(),M}function k(b,v){const{currentType:A}=v;if(A!==7)return!1;j(b);const M=b.currentPeek()===".";return b.resetPeek(),M}function _(b,v){const{currentType:A}=v;if(A!==8)return!1;j(b);const M=g(b.currentPeek());return b.resetPeek(),M}function D(b,v){const{currentType:A}=v;if(!(A===7||A===11))return!1;j(b);const M=b.currentPeek()===":";return b.resetPeek(),M}function S(b,v){const{currentType:A}=v;if(A!==9)return!1;const M=()=>{const K=b.currentPeek();return K==="{"?g(b.peek()):K==="@"||K==="|"||K===":"||K==="."||K===Gn||!K?!1:K===cn?(b.peek(),M()):J(b,!1)},Z=M();return b.resetPeek(),Z}function I(b){j(b);const v=b.currentPeek()==="|";return b.resetPeek(),v}function J(b,v=!0){const A=(Z=!1,K="")=>{const ns=b.currentPeek();return ns==="{"||ns==="@"||!ns?Z:ns==="|"?!(K===Gn||K===cn):ns===Gn?(b.peek(),A(!0,Gn)):ns===cn?(b.peek(),A(!0,cn)):!0},M=A();return v&&b.resetPeek(),M}function U(b,v){const A=b.currentChar();return A===oa?oa:v(A)?(b.next(),A):null}function F(b){const v=b.charCodeAt(0);return v>=97&&v<=122||v>=65&&v<=90||v>=48&&v<=57||v===95||v===36}function G(b){return U(b,F)}function is(b){const v=b.charCodeAt(0);return v>=97&&v<=122||v>=65&&v<=90||v>=48&&v<=57||v===95||v===36||v===45}function ps(b){return U(b,is)}function Y(b){const v=b.charCodeAt(0);return v>=48&&v<=57}function us(b){return U(b,Y)}function Ns(b){const v=b.charCodeAt(0);return v>=48&&v<=57||v>=65&&v<=70||v>=97&&v<=102}function Xs(b){return U(b,Ns)}function Rs(b){let v="",A="";for(;v=us(b);)A+=v;return A}function zs(b){let v="";for(;;){const A=b.currentChar();if(A==="{"||A==="}"||A==="@"||A==="|"||!A)break;if(A===Gn||A===cn)if(J(b))v+=A,b.next();else{if(I(b))break;v+=A,b.next()}else v+=A,b.next()}return v}function ts(b){T(b);let v="",A="";for(;v=ps(b);)A+=v;return b.currentChar()===oa&&i(ms.UNTERMINATED_CLOSING_BRACE,l(),0),A}function rs(b){T(b);let v="";return b.currentChar()==="-"?(b.next(),v+=`-${Rs(b)}`):v+=Rs(b),b.currentChar()===oa&&i(ms.UNTERMINATED_CLOSING_BRACE,l(),0),v}function fs(b){return b!==zp&&b!==cn}function _s(b){T(b),y(b,"'");let v="",A="";for(;v=U(b,fs);)v==="\\"?A+=O(b):A+=v;const M=b.currentChar();return M===cn||M===oa?(i(ms.UNTERMINATED_SINGLE_QUOTE_IN_PLACEHOLDER,l(),0),M===cn&&(b.next(),y(b,"'")),A):(y(b,"'"),A)}function O(b){const v=b.currentChar();switch(v){case"\\":case"'":return b.next(),`\\${v}`;case"u":return W(b,v,4);case"U":return W(b,v,6);default:return i(ms.UNKNOWN_ESCAPE_SEQUENCE,l(),0,v),""}}function W(b,v,A){y(b,v);let M="";for(let Z=0;Z<A;Z++){const K=Xs(b);if(!K){i(ms.INVALID_UNICODE_ESCAPE_SEQUENCE,l(),0,`\\${v}${M}${b.currentChar()}`);break}M+=K}return`\\${v}${M}`}function H(b){return b!=="{"&&b!=="}"&&b!==Gn&&b!==cn}function Q(b){T(b);let v="",A="";for(;v=U(b,H);)A+=v;return A}function gs(b){let v="",A="";for(;v=G(b);)A+=v;return A}function C(b){const v=A=>{const M=b.currentChar();return M==="{"||M==="@"||M==="|"||M==="("||M===")"||!M||M===Gn?A:(A+=M,b.next(),v(A))};return v("")}function P(b){T(b);const v=y(b,"|");return T(b),v}function R(b,v){let A=null;switch(b.currentChar()){case"{":return v.braceNest>=1&&i(ms.NOT_ALLOW_NEST_PLACEHOLDER,l(),0),b.next(),A=d(v,2,"{"),T(b),v.braceNest++,A;case"}":return v.braceNest>0&&v.currentType===2&&i(ms.EMPTY_PLACEHOLDER,l(),0),b.next(),A=d(v,3,"}"),v.braceNest--,v.braceNest>0&&T(b),v.inLinked&&v.braceNest===0&&(v.inLinked=!1),A;case"@":return v.braceNest>0&&i(ms.UNTERMINATED_CLOSING_BRACE,l(),0),A=$(b,v)||h(v),v.braceNest=0,A;default:{let Z=!0,K=!0,ns=!0;if(I(b))return v.braceNest>0&&i(ms.UNTERMINATED_CLOSING_BRACE,l(),0),A=d(v,1,P(b)),v.braceNest=0,v.inLinked=!1,A;if(v.braceNest>0&&(v.currentType===4||v.currentType===5||v.currentType===6))return i(ms.UNTERMINATED_CLOSING_BRACE,l(),0),v.braceNest=0,V(b,v);if(Z=w(b,v))return A=d(v,4,ts(b)),T(b),A;if(K=f(b,v))return A=d(v,5,rs(b)),T(b),A;if(ns=x(b,v))return A=d(v,6,_s(b)),T(b),A;if(!Z&&!K&&!ns)return A=d(v,12,Q(b)),i(ms.INVALID_TOKEN_IN_PLACEHOLDER,l(),0,A.value),T(b),A;break}}return A}function $(b,v){const{currentType:A}=v;let M=null;const Z=b.currentChar();switch((A===7||A===8||A===11||A===9)&&(Z===cn||Z===Gn)&&i(ms.INVALID_LINKED_FORMAT,l(),0),Z){case"@":return b.next(),M=d(v,7,"@"),v.inLinked=!0,M;case".":return T(b),b.next(),d(v,8,".");case":":return T(b),b.next(),d(v,9,":");default:return I(b)?(M=d(v,1,P(b)),v.braceNest=0,v.inLinked=!1,M):k(b,v)||D(b,v)?(T(b),$(b,v)):_(b,v)?(T(b),d(v,11,gs(b))):S(b,v)?(T(b),Z==="{"?R(b,v)||M:d(v,10,C(b))):(A===7&&i(ms.INVALID_LINKED_FORMAT,l(),0),v.braceNest=0,v.inLinked=!1,V(b,v))}}function V(b,v){let A={type:13};if(v.braceNest>0)return R(b,v)||h(v);if(v.inLinked)return $(b,v)||h(v);switch(b.currentChar()){case"{":return R(b,v)||h(v);case"}":return i(ms.UNBALANCED_CLOSING_BRACE,l(),0),b.next(),d(v,3,"}");case"@":return $(b,v)||h(v);default:{if(I(b))return A=d(v,1,P(b)),v.braceNest=0,v.inLinked=!1,A;if(J(b))return d(v,0,zs(b));break}}return A}function q(){const{currentType:b,offset:v,startLoc:A,endLoc:M}=c;return c.lastType=b,c.lastOffset=v,c.lastStartLoc=A,c.lastEndLoc=M,c.offset=t(),c.startLoc=l(),e.currentChar()===oa?d(c,13):V(e,c)}return{nextToken:q,currentOffset:t,currentPosition:l,context:u}}const qj="parser",zj=/(?:\\\\|\\'|\\u([0-9a-fA-F]{4})|\\U([0-9a-fA-F]{6}))/g;function Vj(s,n,a){switch(s){case"\\\\":return"\\";case"\\'":return"'";default:{const e=parseInt(n||a,16);return e<=55295||e>=57344?String.fromCodePoint(e):"�"}}}function Uj(s={}){const n=s.location!==!1,{onError:a}=s;function e(g,m,w,f,...x){const k=g.currentPosition();if(k.offset+=f,k.column+=f,a){const _=n?El(w,k):null,D=Ne(m,_,{domain:qj,args:x});a(D)}}function t(g,m,w){const f={type:g};return n&&(f.start=m,f.end=m,f.loc={start:w,end:w}),f}function l(g,m,w,f){n&&(g.end=m,g.loc&&(g.loc.end=w))}function o(g,m){const w=g.context(),f=t(3,w.offset,w.startLoc);return f.value=m,l(f,g.currentOffset(),g.currentPosition()),f}function p(g,m){const w=g.context(),{lastOffset:f,lastStartLoc:x}=w,k=t(5,f,x);return k.index=parseInt(m,10),g.nextToken(),l(k,g.currentOffset(),g.currentPosition()),k}function c(g,m){const w=g.context(),{lastOffset:f,lastStartLoc:x}=w,k=t(4,f,x);return k.key=m,g.nextToken(),l(k,g.currentOffset(),g.currentPosition()),k}function u(g,m){const w=g.context(),{lastOffset:f,lastStartLoc:x}=w,k=t(9,f,x);return k.value=m.replace(zj,Vj),g.nextToken(),l(k,g.currentOffset(),g.currentPosition()),k}function r(g){const m=g.nextToken(),w=g.context(),{lastOffset:f,lastStartLoc:x}=w,k=t(8,f,x);return m.type!==11?(e(g,ms.UNEXPECTED_EMPTY_LINKED_MODIFIER,w.lastStartLoc,0),k.value="",l(k,f,x),{nextConsumeToken:m,node:k}):(m.value==null&&e(g,ms.UNEXPECTED_LEXICAL_ANALYSIS,w.lastStartLoc,0,In(m)),k.value=m.value||"",l(k,g.currentOffset(),g.currentPosition()),{node:k})}function i(g,m){const w=g.context(),f=t(7,w.offset,w.startLoc);return f.value=m,l(f,g.currentOffset(),g.currentPosition()),f}function d(g){const m=g.context(),w=t(6,m.offset,m.startLoc);let f=g.nextToken();if(f.type===8){const x=r(g);w.modifier=x.node,f=x.nextConsumeToken||g.nextToken()}switch(f.type!==9&&e(g,ms.UNEXPECTED_LEXICAL_ANALYSIS,m.lastStartLoc,0,In(f)),f=g.nextToken(),f.type===2&&(f=g.nextToken()),f.type){case 10:f.value==null&&e(g,ms.UNEXPECTED_LEXICAL_ANALYSIS,m.lastStartLoc,0,In(f)),w.key=i(g,f.value||"");break;case 4:f.value==null&&e(g,ms.UNEXPECTED_LEXICAL_ANALYSIS,m.lastStartLoc,0,In(f)),w.key=c(g,f.value||"");break;case 5:f.value==null&&e(g,ms.UNEXPECTED_LEXICAL_ANALYSIS,m.lastStartLoc,0,In(f)),w.key=p(g,f.value||"");break;case 6:f.value==null&&e(g,ms.UNEXPECTED_LEXICAL_ANALYSIS,m.lastStartLoc,0,In(f)),w.key=u(g,f.value||"");break;default:{e(g,ms.UNEXPECTED_EMPTY_LINKED_KEY,m.lastStartLoc,0);const x=g.context(),k=t(7,x.offset,x.startLoc);return k.value="",l(k,x.offset,x.startLoc),w.key=k,l(w,x.offset,x.startLoc),{nextConsumeToken:f,node:w}}}return l(w,g.currentOffset(),g.currentPosition()),{node:w}}function h(g){const m=g.context(),w=m.currentType===1?g.currentOffset():m.offset,f=m.currentType===1?m.endLoc:m.startLoc,x=t(2,w,f);x.items=[];let k=null;do{const S=k||g.nextToken();switch(k=null,S.type){case 0:S.value==null&&e(g,ms.UNEXPECTED_LEXICAL_ANALYSIS,m.lastStartLoc,0,In(S)),x.items.push(o(g,S.value||""));break;case 5:S.value==null&&e(g,ms.UNEXPECTED_LEXICAL_ANALYSIS,m.lastStartLoc,0,In(S)),x.items.push(p(g,S.value||""));break;case 4:S.value==null&&e(g,ms.UNEXPECTED_LEXICAL_ANALYSIS,m.lastStartLoc,0,In(S)),x.items.push(c(g,S.value||""));break;case 6:S.value==null&&e(g,ms.UNEXPECTED_LEXICAL_ANALYSIS,m.lastStartLoc,0,In(S)),x.items.push(u(g,S.value||""));break;case 7:{const I=d(g);x.items.push(I.node),k=I.nextConsumeToken||null;break}}}while(m.currentType!==13&&m.currentType!==1);const _=m.currentType===1?m.lastOffset:g.currentOffset(),D=m.currentType===1?m.lastEndLoc:g.currentPosition();return l(x,_,D),x}function y(g,m,w,f){const x=g.context();let k=f.items.length===0;const _=t(1,m,w);_.cases=[],_.cases.push(f);do{const D=h(g);k||(k=D.items.length===0),_.cases.push(D)}while(x.currentType!==13);return k&&e(g,ms.MUST_HAVE_MESSAGES_IN_PLURAL,w,0),l(_,g.currentOffset(),g.currentPosition()),_}function j(g){const m=g.context(),{offset:w,startLoc:f}=m,x=h(g);return m.currentType===13?x:y(g,w,f,x)}function T(g){const m=Bj(g,sn({},s)),w=m.context(),f=t(0,w.offset,w.startLoc);return n&&f.loc&&(f.loc.source=g),f.body=j(m),s.onCacheKey&&(f.cacheKey=s.onCacheKey(g)),w.currentType!==13&&e(m,ms.UNEXPECTED_LEXICAL_ANALYSIS,w.lastStartLoc,0,g[w.offset]||""),l(f,m.currentOffset(),m.currentPosition()),f}return{parse:T}}function In(s){if(s.type===13)return"EOF";const n=(s.value||"").replace(/\r?\n/gu,"\\n");return n.length>10?n.slice(0,9)+"…":n}const Hj=/<\/?[\w\s="/.':;#-\/]+>/,Gj=s=>Hj.test(s);function En(s){return ks(s)&&fo(s)===0&&(An(s,"b")||An(s,"body"))}const ci=["b","body"];function Wj(s){return ya(s,ci)}const ri=["c","cases"];function Kj(s){return ya(s,ri,[])}const ii=["s","static"];function Xj(s){return ya(s,ii)}const ui=["i","items"];function Yj(s){return ya(s,ui,[])}const hi=["t","type"];function fo(s){return ya(s,hi)}const di=["v","value"];function We(s,n){const a=ya(s,di);if(a!=null)return a;throw ke(n)}const mi=["m","modifier"];function Qj(s){return ya(s,mi)}const ji=["k","key"];function Jj(s){const n=ya(s,ji);if(n)return n;throw ke(6)}function ya(s,n,a){for(let e=0;e<n.length;e++){const t=n[e];if(An(s,t)&&s[t]!=null)return s[t]}return a}const fi=[...ci,...ri,...ii,...ui,...ji,...mi,...di,...hi];function ke(s){return new Error(`unhandled node type: ${s}`)}function al(s){return a=>Zj(a,s)}function Zj(s,n){const a=Wj(n);if(a==null)throw ke(0);if(fo(a)===1){const l=Kj(a);return s.plural(l.reduce((o,p)=>[...o,Vp(s,p)],[]))}else return Vp(s,a)}function Vp(s,n){const a=Xj(n);if(a!=null)return s.type==="text"?a:s.normalize([a]);{const e=Yj(n).reduce((t,l)=>[...t,Sl(s,l)],[]);return s.normalize(e)}}function Sl(s,n){const a=fo(n);switch(a){case 3:return We(n,a);case 9:return We(n,a);case 4:{const e=n;if(An(e,"k")&&e.k)return s.interpolate(s.named(e.k));if(An(e,"key")&&e.key)return s.interpolate(s.named(e.key));throw ke(a)}case 5:{const e=n;if(An(e,"i")&&Ws(e.i))return s.interpolate(s.list(e.i));if(An(e,"index")&&Ws(e.index))return s.interpolate(s.list(e.index));throw ke(a)}case 6:{const e=n,t=Qj(e),l=Jj(e);return s.linked(Sl(s,l),t?Sl(s,t):void 0,s.type)}case 7:return We(n,a);case 8:return We(n,a);default:throw new Error(`unhandled node on format message part: ${a}`)}}const sf="Detected HTML in '{source}' message. Recommend not using HTML messages to avoid XSS.";function nf(s,n){n&&Gj(s)&&_a(Mt(sf,{source:s}))}const af=s=>s;let Ke=Ts();function ef(s,n={}){let a=!1;const e=n.onError||Rj;n.onError=o=>{a=!0,e(o)};const l=Uj(n).parse(s);return n.optimize&&Mj(l),n.mangle&&$a(l),{ast:l,detectError:a,code:""}}function tf(s,n){if(es(s)){const a=Is(n.warnHtmlMessage)?n.warnHtmlMessage:!0;nf(s,a);const t=(n.onCacheKey||af)(s),l=Ke[t];if(l)return l;const{ast:o,detectError:p}=ef(s,{...n,location:!0,mangle:!1,optimize:!1}),c=al(o);return p?c:Ke[t]=c}else{if(!En(s))return _a(`the message that is resolve with key '${n.key}' is not supported for jit compilation`),(()=>s);const a=s.cacheKey;if(a){const e=Ke[a];return e||(Ke[a]=al(s))}else return al(s)}}let xe=null;function lf(s){xe=s}function of(s,n,a){xe&&xe.emit("i18n:init",{timestamp:Date.now(),i18n:s,version:n,meta:a})}const pf=cf("function:translate");function cf(s){return n=>xe&&xe.emit(s,n)}const tn={INVALID_ARGUMENT:Tj,INVALID_DATE_ARGUMENT:18,INVALID_ISO_DATE_ARGUMENT:19,NOT_SUPPORT_NON_STRING_MESSAGE:20,NOT_SUPPORT_LOCALE_PROMISE_VALUE:21,NOT_SUPPORT_LOCALE_ASYNC_FUNCTION:22,NOT_SUPPORT_LOCALE_TYPE:23},rf=24;function Jn(s){return Ne(s,null,{messages:uf})}const uf={[tn.INVALID_ARGUMENT]:"Invalid arguments",[tn.INVALID_DATE_ARGUMENT]:"The date provided is an invalid Date object.Make sure your Date represents a valid date.",[tn.INVALID_ISO_DATE_ARGUMENT]:"The argument provided is not a valid ISO date string",[tn.NOT_SUPPORT_NON_STRING_MESSAGE]:"Not support non-string message",[tn.NOT_SUPPORT_LOCALE_PROMISE_VALUE]:"cannot support promise value",[tn.NOT_SUPPORT_LOCALE_ASYNC_FUNCTION]:"cannot support async function",[tn.NOT_SUPPORT_LOCALE_TYPE]:"cannot support locale type"};function go(s,n){return n.locale!=null?Up(n.locale):Up(s.locale)}let el;function Up(s){if(es(s))return s;if(Os(s)){if(s.resolvedOnce&&el!=null)return el;if(s.constructor.name==="Function"){const n=s();if(kj(n))throw Jn(tn.NOT_SUPPORT_LOCALE_PROMISE_VALUE);return el=n}else throw Jn(tn.NOT_SUPPORT_LOCALE_ASYNC_FUNCTION)}else throw Jn(tn.NOT_SUPPORT_LOCALE_TYPE)}function hf(s,n,a){return[...new Set([a,...qs(n)?n:ks(n)?Object.keys(n):es(n)?[n]:[a]])]}function gi(s,n,a){const e=es(a)?a:gt,t=s;t.__localeChainCache||(t.__localeChainCache=new Map);let l=t.__localeChainCache.get(e);if(!l){l=[];let o=[a];for(;qs(o);)o=Hp(l,o,n);const p=qs(n)||!Cs(n)?n:n.default?n.default:null;o=es(p)?[p]:p,qs(o)&&Hp(l,o,!1),t.__localeChainCache.set(e,l)}return l}function Hp(s,n,a){let e=!0;for(let t=0;t<n.length&&Is(e);t++){const l=n[t];es(l)&&(e=df(s,n[t],a))}return e}function df(s,n,a){let e;const t=n.split("-");do{const l=t.join("-");e=mf(s,l,a),t.splice(-1,1)}while(t.length&&e===!0);return e}function mf(s,n,a){let e=!1;if(!s.includes(n)&&(e=!0,n)){e=n[n.length-1]!=="!";const t=n.replace(/!/g,"");s.push(t),(qs(a)||Cs(a))&&a[t]&&(e=a[t])}return e}const va=[];va[0]={w:[0],i:[3,0],"[":[4],o:[7]};va[1]={w:[1],".":[2],"[":[4],o:[7]};va[2]={w:[2],i:[3,0],0:[3,0]};va[3]={i:[3,0],0:[3,0],w:[1,1],".":[2,1],"[":[4,1],o:[7,1]};va[4]={"'":[5,0],'"':[6,0],"[":[4,2],"]":[1,3],o:8,l:[4,0]};va[5]={"'":[4,0],o:8,l:[5,0]};va[6]={'"':[4,0],o:8,l:[6,0]};const jf=/^\s?(?:true|false|-?[\d.]+|'[^']*'|"[^"]*")\s?$/;function ff(s){return jf.test(s)}function gf(s){const n=s.charCodeAt(0),a=s.charCodeAt(s.length-1);return n===a&&(n===34||n===39)?s.slice(1,-1):s}function bf(s){if(s==null)return"o";switch(s.charCodeAt(0)){case 91:case 93:case 46:case 34:case 39:return s;case 95:case 36:case 45:return"i";case 9:case 10:case 13:case 160:case 65279:case 8232:case 8233:return"w"}return"i"}function _f(s){const n=s.trim();return s.charAt(0)==="0"&&isNaN(parseInt(s))?!1:ff(n)?gf(n):"*"+n}function yf(s){const n=[];let a=-1,e=0,t=0,l,o,p,c,u,r,i;const d=[];d[0]=()=>{o===void 0?o=p:o+=p},d[1]=()=>{o!==void 0&&(n.push(o),o=void 0)},d[2]=()=>{d[0](),t++},d[3]=()=>{if(t>0)t--,e=4,d[0]();else{if(t=0,o===void 0||(o=_f(o),o===!1))return!1;d[1]()}};function h(){const y=s[a+1];if(e===5&&y==="'"||e===6&&y==='"')return a++,p="\\"+y,d[0](),!0}for(;e!==null;)if(a++,l=s[a],!(l==="\\"&&h())){if(c=bf(l),i=va[e],u=i[c]||i.l||8,u===8||(e=u[0],u[1]!==void 0&&(r=d[u[1]],r&&(p=l,r()===!1))))return;if(e===7)return n}}const Gp=new Map;function vf(s,n){return ks(s)?s[n]:null}function wf(s,n){if(!ks(s))return null;let a=Gp.get(n);if(a||(a=yf(n),a&&Gp.set(n,a)),!a)return null;const e=a.length;let t=s,l=0;for(;l<e;){const o=a[l];if(fi.includes(o)&&En(t))return null;const p=t[o];if(p===void 0||Os(t))return null;t=p,l++}return t}const mn={NOT_FOUND_KEY:1,FALLBACK_TO_TRANSLATE:2,CANNOT_FORMAT_NUMBER:3,FALLBACK_TO_NUMBER_FORMAT:4,CANNOT_FORMAT_DATE:5,FALLBACK_TO_DATE_FORMAT:6,EXPERIMENTAL_CUSTOM_MESSAGE_COMPILER:7},Cf=8,kf={[mn.NOT_FOUND_KEY]:"Not found '{key}' key in '{locale}' locale messages.",[mn.FALLBACK_TO_TRANSLATE]:"Fall back to translate '{key}' key with '{target}' locale.",[mn.CANNOT_FORMAT_NUMBER]:"Cannot format a number value due to not supported Intl.NumberFormat.",[mn.FALLBACK_TO_NUMBER_FORMAT]:"Fall back to number format '{key}' key with '{target}' locale.",[mn.CANNOT_FORMAT_DATE]:"Cannot format a date value due to not supported Intl.DateTimeFormat.",[mn.FALLBACK_TO_DATE_FORMAT]:"Fall back to datetime format '{key}' key with '{target}' locale.",[mn.EXPERIMENTAL_CUSTOM_MESSAGE_COMPILER]:"This project is using Custom Message Compiler, which is an experimental feature. It may receive breaking changes or be removed in the future."};function La(s,...n){return Mt(kf[s],...n)}const xf="12.0.0-alpha.3",Ot=-1,gt="en-US",bt="",Wp=s=>`${s.charAt(0).toLocaleUpperCase()}${s.substr(1)}`;function Pf(){return{upper:(s,n)=>n==="text"&&es(s)?s.toUpperCase():n==="vnode"&&ks(s)&&"__v_isVNode"in s?s.children.toUpperCase():s,lower:(s,n)=>n==="text"&&es(s)?s.toLowerCase():n==="vnode"&&ks(s)&&"__v_isVNode"in s?s.children.toLowerCase():s,capitalize:(s,n)=>n==="text"&&es(s)?Wp(s):n==="vnode"&&ks(s)&&"__v_isVNode"in s?Wp(s.children):s}}let bi;function Ef(s){bi=s}let _i;function Sf(s){_i=s}let yi;function Tf(s){yi=s}let vi=null;const Af=s=>{vi=s},Rf=()=>vi;let wi=null;const Kp=s=>{wi=s},Lf=()=>wi;let Xp=0;function Mf(s={}){const n=Os(s.onWarn)?s.onWarn:_a,a=es(s.version)?s.version:xf,e=es(s.locale)||Os(s.locale)?s.locale:gt,t=Os(e)?gt:e,l=qs(s.fallbackLocale)||Cs(s.fallbackLocale)||es(s.fallbackLocale)||s.fallbackLocale===!1?s.fallbackLocale:t,o=Cs(s.messages)?s.messages:tl(t),p=Cs(s.datetimeFormats)?s.datetimeFormats:tl(t),c=Cs(s.numberFormats)?s.numberFormats:tl(t),u=sn(Ts(),s.modifiers,Pf()),r=s.pluralRules||Ts(),i=Os(s.missing)?s.missing:null,d=Is(s.missingWarn)||ft(s.missingWarn)?s.missingWarn:!0,h=Is(s.fallbackWarn)||ft(s.fallbackWarn)?s.fallbackWarn:!0,y=!!s.fallbackFormat,j=!!s.unresolving,T=Os(s.postTranslation)?s.postTranslation:null,g=Cs(s.processor)?s.processor:null,m=Is(s.warnHtmlMessage)?s.warnHtmlMessage:!0,w=!!s.escapeParameter,f=Os(s.messageCompiler)?s.messageCompiler:bi;Os(s.messageCompiler)&&Ej(La(mn.EXPERIMENTAL_CUSTOM_MESSAGE_COMPILER));const x=Os(s.messageResolver)?s.messageResolver:_i||vf,k=Os(s.localeFallbacker)?s.localeFallbacker:yi||hf,_=ks(s.fallbackContext)?s.fallbackContext:void 0,D=s,S=ks(D.__datetimeFormatters)?D.__datetimeFormatters:new Map,I=ks(D.__numberFormatters)?D.__numberFormatters:new Map,J=ks(D.__meta)?D.__meta:{};Xp++;const U={version:a,cid:Xp,locale:e,fallbackLocale:l,messages:o,modifiers:u,pluralRules:r,missing:i,missingWarn:d,fallbackWarn:h,fallbackFormat:y,unresolving:j,postTranslation:T,processor:g,warnHtmlMessage:m,escapeParameter:w,messageCompiler:f,messageResolver:x,localeFallbacker:k,fallbackContext:_,onWarn:n,__meta:J};return U.datetimeFormats=p,U.numberFormats=c,U.__datetimeFormatters=S,U.__numberFormatters=I,U.__v_emitter=D.__v_emitter!=null?D.__v_emitter:void 0,of(U,a,J),U}const tl=s=>({[s]:Ts()});function It(s,n){return s instanceof RegExp?s.test(n):s}function Ci(s,n){return s instanceof RegExp?s.test(n):s}function bo(s,n,a,e,t){const{missing:l,onWarn:o}=s;{const p=s.__v_emitter;p&&p.emit("missing",{locale:a,key:n,type:t,groupId:`${t}:${n}`})}if(l!==null){const p=l(s,a,n,t);return es(p)?p:n}else return Ci(e,n)&&o(La(mn.NOT_FOUND_KEY,{key:n,locale:a})),n}function ne(s,n,a){const e=s;e.__localeChainCache=new Map,s.localeFallbacker(s,a,n)}function ki(s,n){return s===n?!1:s.split("-")[0]===n.split("-")[0]}function Df(s,n){const a=n.indexOf(s);if(a===-1)return!1;for(let e=a+1;e<n.length;e++)if(ki(s,n[e]))return!0;return!1}const Yp=typeof Intl<"u",xi={dateTimeFormat:Yp&&typeof Intl.DateTimeFormat<"u",numberFormat:Yp&&typeof Intl.NumberFormat<"u"};function Qp(s,...n){const{datetimeFormats:a,unresolving:e,fallbackLocale:t,onWarn:l,localeFallbacker:o}=s,{__datetimeFormatters:p}=s;if(!xi.dateTimeFormat)return l(La(mn.CANNOT_FORMAT_DATE)),bt;const[c,u,r,i]=Tl(...n),d=Is(r.missingWarn)?r.missingWarn:s.missingWarn,h=Is(r.fallbackWarn)?r.fallbackWarn:s.fallbackWarn,y=!!r.part,j=go(s,r),T=o(s,t,j);if(!es(c)||c==="")return new Intl.DateTimeFormat(j,i).format(u);let g={},m,w=null,f=j,x=null;const k="datetime format";for(let S=0;S<T.length;S++){if(m=x=T[S],j!==m&&It(h,c)&&l(La(mn.FALLBACK_TO_DATE_FORMAT,{key:c,target:m})),j!==m){const I=s.__v_emitter;I&&I.emit("fallback",{type:k,key:c,from:f,to:x,groupId:`${k}:${c}`})}if(g=a[m]||{},w=g[c],Cs(w))break;bo(s,c,m,d,k),f=x}if(!Cs(w)||!es(m))return e?Ot:c;let _=`${m}__${c}`;Dt(i)||(_=`${_}__${JSON.stringify(i)}`);let D=p.get(_);return D||(D=new Intl.DateTimeFormat(m,sn({},w,i)),p.set(_,D)),y?D.formatToParts(u):D.format(u)}const Pi=["localeMatcher","weekday","era","year","month","day","hour","minute","second","timeZoneName","formatMatcher","hour12","timeZone","dateStyle","timeStyle","calendar","dayPeriod","numberingSystem","hourCycle","fractionalSecondDigits"];function Tl(...s){const[n,a,e,t]=s,l=Ts();let o=Ts(),p;if(es(n)){const c=n.match(/(\d{4}-\d{2}-\d{2})(T|\s)?(.*)/);if(!c)throw Jn(tn.INVALID_ISO_DATE_ARGUMENT);const u=c[3]?c[3].trim().startsWith("T")?`${c[1].trim()}${c[3].trim()}`:`${c[1].trim()}T${c[3].trim()}`:c[1].trim();p=new Date(u);try{p.toISOString()}catch{throw Jn(tn.INVALID_ISO_DATE_ARGUMENT)}}else if(yj(n)){if(isNaN(n.getTime()))throw Jn(tn.INVALID_DATE_ARGUMENT);p=n}else if(Ws(n))p=n;else throw Jn(tn.INVALID_ARGUMENT);return es(a)?l.key=a:Cs(a)&&Object.keys(a).forEach(c=>{Pi.includes(c)?o[c]=a[c]:l[c]=a[c]}),es(e)?l.locale=e:Cs(e)&&(o=e),Cs(t)&&(o=t),[l.key||"",p,l,o]}function Jp(s,n,a){const e=s;for(const t in a){const l=`${n}__${t}`;e.__datetimeFormatters.has(l)&&e.__datetimeFormatters.delete(l)}}function Zp(s,...n){const{numberFormats:a,unresolving:e,fallbackLocale:t,onWarn:l,localeFallbacker:o}=s,{__numberFormatters:p}=s;if(!xi.numberFormat)return l(La(mn.CANNOT_FORMAT_NUMBER)),bt;const[c,u,r,i]=Al(...n),d=Is(r.missingWarn)?r.missingWarn:s.missingWarn,h=Is(r.fallbackWarn)?r.fallbackWarn:s.fallbackWarn,y=!!r.part,j=go(s,r),T=o(s,t,j);if(!es(c)||c==="")return new Intl.NumberFormat(j,i).format(u);let g={},m,w=null,f=j,x=null;const k="number format";for(let S=0;S<T.length;S++){if(m=x=T[S],j!==m&&It(h,c)&&l(La(mn.FALLBACK_TO_NUMBER_FORMAT,{key:c,target:m})),j!==m){const I=s.__v_emitter;I&&I.emit("fallback",{type:k,key:c,from:f,to:x,groupId:`${k}:${c}`})}if(g=a[m]||{},w=g[c],Cs(w))break;bo(s,c,m,d,k),f=x}if(!Cs(w)||!es(m))return e?Ot:c;let _=`${m}__${c}`;Dt(i)||(_=`${_}__${JSON.stringify(i)}`);let D=p.get(_);return D||(D=new Intl.NumberFormat(m,sn({},w,i)),p.set(_,D)),y?D.formatToParts(u):D.format(u)}const Ei=["localeMatcher","style","currency","currencyDisplay","currencySign","useGrouping","minimumIntegerDigits","minimumFractionDigits","maximumFractionDigits","minimumSignificantDigits","maximumSignificantDigits","compactDisplay","notation","signDisplay","unit","unitDisplay","roundingMode","roundingPriority","roundingIncrement","trailingZeroDisplay"];function Al(...s){const[n,a,e,t]=s,l=Ts();let o=Ts();if(!Ws(n))throw Jn(tn.INVALID_ARGUMENT);const p=n;return es(a)?l.key=a:Cs(a)&&Object.keys(a).forEach(c=>{Ei.includes(c)?o[c]=a[c]:l[c]=a[c]}),es(e)?l.locale=e:Cs(e)&&(o=e),Cs(t)&&(o=t),[l.key||"",p,l,o]}function sc(s,n,a){const e=s;for(const t in a){const l=`${n}__${t}`;e.__numberFormatters.has(l)&&e.__numberFormatters.delete(l)}}const Of=s=>s,If=s=>"",Nf="text",Ff=s=>s.length===0?"":oi(s),$f=xj;function nc(s,n){return s=Math.abs(s),n===2?s?s>1?1:0:1:s?Math.min(s,2):0}function Bf(s){const n=Ws(s.pluralIndex)?s.pluralIndex:-1;return s.named&&(Ws(s.named.count)||Ws(s.named.n))?Ws(s.named.count)?s.named.count:Ws(s.named.n)?s.named.n:n:n}function qf(s,n){n.count||(n.count=s),n.n||(n.n=s)}function zf(s={}){const n=s.locale,a=Bf(s),e=ks(s.pluralRules)&&es(n)&&Os(s.pluralRules[n])?s.pluralRules[n]:nc,t=ks(s.pluralRules)&&es(n)&&Os(s.pluralRules[n])?nc:void 0,l=g=>g[e(a,g.length,t)],o=s.list||[],p=g=>o[g],c=s.named||Ts();Ws(s.pluralIndex)&&qf(a,c);const u=g=>c[g];function r(g,m){const w=Os(s.messages)?s.messages(g,!!m):ks(s.messages)?s.messages[g]:!1;return w||(s.parent?s.parent.message(g):If)}const i=g=>s.modifiers?s.modifiers[g]:Of,d=Cs(s.processor)&&Os(s.processor.normalize)?s.processor.normalize:Ff,h=Cs(s.processor)&&Os(s.processor.interpolate)?s.processor.interpolate:$f,y=Cs(s.processor)&&es(s.processor.type)?s.processor.type:Nf,T={list:p,named:u,plural:l,linked:(g,...m)=>{const[w,f]=m;let x="text",k="";m.length===1?ks(w)?(k=w.modifier||k,x=w.type||x):es(w)&&(k=w||k):m.length===2&&(es(w)&&(k=w||k),es(f)&&(x=f||x));const _=r(g,!0)(T),D=x==="vnode"&&qs(_)&&k?_[0]:_;return k?i(k)(D,x):D},message:r,type:y,interpolate:h,normalize:d,values:sn(Ts(),o,c)};return T}const ac=()=>"",xn=s=>Os(s);function ec(s,...n){const{fallbackFormat:a,postTranslation:e,unresolving:t,messageCompiler:l,fallbackLocale:o,messages:p}=s,[c,u]=Rl(...n),r=Is(u.missingWarn)?u.missingWarn:s.missingWarn,i=Is(u.fallbackWarn)?u.fallbackWarn:s.fallbackWarn,d=Is(u.escapeParameter)?u.escapeParameter:s.escapeParameter,h=!!u.resolvedMessage,y=es(u.default)||Is(u.default)?Is(u.default)?l?c:()=>c:u.default:a?l?c:()=>c:null,j=a||y!=null&&(es(y)||Os(y)),T=go(s,u);d&&Vf(u);let[g,m,w]=h?[c,T,p[T]||Ts()]:Si(s,c,T,o,i,r),f=g,x=c;if(!h&&!(es(f)||En(f)||xn(f))&&j&&(f=y,x=f),!h&&(!(es(f)||En(f)||xn(f))||!es(m)))return t?Ot:c;if(es(f)&&s.messageCompiler==null)return _a(`The message format compilation is not supported in this build. Because message compiler isn't included. You need to pre-compilation all message format. So translate function return '${c}'.`),c;let k=!1;const _=()=>{k=!0},D=xn(f)?f:Ti(s,c,m,f,x,_);if(k)return f;const S=Wf(s,m,w,u),I=zf(S),J=Uf(s,D,I),U=e?e(J,c):J;{const F={timestamp:Date.now(),key:es(c)?c:xn(f)?f.key:"",locale:m||(xn(f)?f.locale:""),format:es(f)?f:xn(f)?f.source:"",message:U};F.meta=sn({},s.__meta,Rf()||{}),pf(F)}return U}function Vf(s){qs(s.list)?s.list=s.list.map(n=>es(n)?Fp(n):n):ks(s.named)&&Object.keys(s.named).forEach(n=>{es(s.named[n])&&(s.named[n]=Fp(s.named[n]))})}function Si(s,n,a,e,t,l){const{messages:o,onWarn:p,messageResolver:c,localeFallbacker:u}=s,r=u(s,e,a);let i=Ts(),d,h=null,y=a,j=null;const T="translate";for(let g=0;g<r.length;g++){if(d=j=r[g],a!==d&&!ki(a,d)&&It(t,n)&&p(La(mn.FALLBACK_TO_TRANSLATE,{key:n,target:d})),a!==d){const x=s.__v_emitter;x&&x.emit("fallback",{type:T,key:n,from:y,to:j,groupId:`${T}:${n}`})}i=o[d]||Ts();let m=null,w,f;if(na&&(m=window.performance.now(),w="intlify-message-resolve-start",f="intlify-message-resolve-end",yn&&yn(w)),(h=c(i,n))===null&&(h=i[n]),na){const x=window.performance.now(),k=s.__v_emitter;k&&m&&h&&k.emit("message-resolve",{type:"message-resolve",key:n,message:h,time:x-m,groupId:`${T}:${n}`}),w&&f&&yn&&Ra&&(yn(f),Ra("intlify message resolve",w,f))}if(es(h)||En(h)||xn(h))break;if(!Df(d,r)){const x=bo(s,n,d,l,T);x!==n&&(h=x)}y=j}return[h,d,i]}function Ti(s,n,a,e,t,l){const{messageCompiler:o,warnHtmlMessage:p}=s;if(xn(e)){const d=e;return d.locale=d.locale||a,d.key=d.key||n,d}if(o==null){const d=(()=>e);return d.locale=a,d.key=n,d}let c=null,u,r;na&&(c=window.performance.now(),u="intlify-message-compilation-start",r="intlify-message-compilation-end",yn&&yn(u));const i=o(e,Hf(s,a,t,e,p,l));if(na){const d=window.performance.now(),h=s.__v_emitter;h&&c&&h.emit("message-compilation",{type:"message-compilation",message:e,time:d-c,groupId:`translate:${n}`}),u&&r&&yn&&Ra&&(yn(r),Ra("intlify message compilation",u,r))}return i.locale=a,i.key=n,i.source=e,i}function Uf(s,n,a){let e=null,t,l;na&&(e=window.performance.now(),t="intlify-message-evaluation-start",l="intlify-message-evaluation-end",yn&&yn(t));const o=n(a);if(na){const p=window.performance.now(),c=s.__v_emitter;c&&e&&c.emit("message-evaluation",{type:"message-evaluation",value:o,time:p-e,groupId:`translate:${n.key}`}),t&&l&&yn&&Ra&&(yn(l),Ra("intlify message evaluation",t,l))}return o}function Rl(...s){const[n,a,e]=s,t=Ts();if(!es(n)&&!Ws(n)&&!xn(n)&&!En(n))throw Jn(tn.INVALID_ARGUMENT);const l=Ws(n)?String(n):(xn(n),n);return Ws(a)?t.plural=a:es(a)?t.default=a:Cs(a)&&!Dt(a)?t.named=a:qs(a)&&(t.list=a),Ws(e)?t.plural=e:es(e)?t.default=e:Cs(e)&&sn(t,e),[l,t]}function Hf(s,n,a,e,t,l){return{locale:n,key:a,warnHtmlMessage:t,onError:o=>{l&&l(o);{const p=Gf(e),c=`Message compilation error: ${o.message}`,u=o.location&&p&&Pj(p,o.location.start.offset,o.location.end.offset),r=s.__v_emitter;r&&p&&r.emit("compile-error",{message:p,error:o.message,start:o.location&&o.location.start.offset,end:o.location&&o.location.end.offset,groupId:`translate:${a}`}),console.error(u?`${c}
${u}`:c)}},onCacheKey:o=>bj(n,a,o)}}function Gf(s){if(es(s))return s;if(s.loc&&s.loc.source)return s.loc.source}function Wf(s,n,a,e){const{modifiers:t,pluralRules:l,messageResolver:o,fallbackLocale:p,fallbackWarn:c,missingWarn:u,fallbackContext:r}=s,d={locale:n,modifiers:t,pluralRules:l,messages:(h,y)=>{let j=o(a,h);if(j==null&&(r||y)){const[,,T]=Si(r||s,h,n,p,c,u);j=o(T,h)}if(es(j)||En(j)){let T=!1;const m=Ti(s,h,n,j,h,()=>{T=!0});return T?ac:m}else return xn(j)?j:ac}};return s.processor&&(d.processor=s.processor),e.list&&(d.list=e.list),e.named&&(d.named=e.named),Ws(e.plural)&&(d.pluralIndex=e.plural),d}const Kf="12.0.0-alpha.3";function Xf(){console.info(`You are running a development build of vue-i18n.
Make sure to use the production build (*.prod.js) when deploying for production.`)}const Vs={UNEXPECTED_RETURN_TYPE:rf,INVALID_ARGUMENT:25,MUST_BE_CALL_SETUP_TOP:26,NOT_INSTALLED:27,REQUIRED_VALUE:28,INVALID_VALUE:29,CANNOT_SETUP_VUE_DEVTOOLS_PLUGIN:30,NOT_INSTALLED_WITH_PROVIDE:31,UNEXPECTED_ERROR:32,DUPLICATE_USE_I18N_CALLING:33};function Bn(s,...n){return Ne(s,null,{messages:Yf,args:n})}const Yf={[Vs.UNEXPECTED_RETURN_TYPE]:"Unexpected return type in composer",[Vs.INVALID_ARGUMENT]:"Invalid argument",[Vs.MUST_BE_CALL_SETUP_TOP]:"Must be called at the top of a `setup` function",[Vs.NOT_INSTALLED]:"Need to install with `app.use` function",[Vs.UNEXPECTED_ERROR]:"Unexpected error",[Vs.REQUIRED_VALUE]:"Required in value: {0}",[Vs.INVALID_VALUE]:"Invalid value",[Vs.CANNOT_SETUP_VUE_DEVTOOLS_PLUGIN]:"Cannot setup vue-devtools plugin",[Vs.NOT_INSTALLED_WITH_PROVIDE]:"Need to install with `provide` function",[Vs.DUPLICATE_USE_I18N_CALLING]:"Duplicate local-scope `useI18n` call detected. Call `useI18n` only once per component."},Ll=zn("__translateVNode"),Ml=zn("__datetimeParts"),Dl=zn("__numberParts"),Pe=zn("__enableEmitter"),Ol=zn("__disableEmitter"),Qf=zn("__setPluralRules"),Jf=zn("__injectWithOption"),Il=zn("__dispose"),Va={FALLBACK_TO_ROOT:Cf,NOT_FOUND_PARENT_SCOPE:9,IGNORE_OBJ_FLATTEN:10},Zf={[Va.FALLBACK_TO_ROOT]:"Fall back to {type} '{key}' with root locale.",[Va.NOT_FOUND_PARENT_SCOPE]:"Not found parent scope. use the global scope.",[Va.IGNORE_OBJ_FLATTEN]:"Ignore object flatten: '{key}' key has an string value"};function _o(s,...n){return Mt(Zf[s],...n)}function Ee(s){if(!ks(s)||En(s))return s;for(const n in s)if(An(s,n))if(!n.includes("."))ks(s[n])&&Ee(s[n]);else{const a=n.split("."),e=a.length-1;let t=s,l=!1;for(let o=0;o<e;o++){if(a[o]==="__proto__")throw new Error(`unsafe key: ${a[o]}`);if(a[o]in t||(t[a[o]]=Ts()),!ks(t[a[o]])){_a(_o(Va.IGNORE_OBJ_FLATTEN,{key:a[o]})),l=!0;break}t=t[a[o]]}if(l||(En(t)?fi.includes(a[e])||delete s[n]:(t[a[e]]=s[n],delete s[n])),!En(t)){const o=t[a[e]];ks(o)&&Ee(o)}}return s}function Ai(s,n){const{messages:a,__i18n:e,messageResolver:t,flatJson:l}=n,o=Cs(a)?a:qs(e)?Ts():{[s]:Ts()};if(qs(e)&&e.forEach(p=>{if("locale"in p&&"resource"in p){const{locale:c,resource:u}=p;c?(o[c]=o[c]||Ts(),et(u,o[c])):et(u,o)}else es(p)&&et(JSON.parse(p),o)}),t==null&&l)for(const p in o)An(o,p)&&Ee(o[p]);return o}function Ri(s){return s.type}function sg(s,n,a){let e=ks(n.messages)?n.messages:Ts();"__i18nGlobal"in a&&(e=Ai(s.locale.value,{messages:e,__i18n:a.__i18nGlobal}));const t=Object.keys(e);t.length&&t.forEach(l=>{s.mergeLocaleMessage(l,e[l])});{if(ks(n.datetimeFormats)){const l=Object.keys(n.datetimeFormats);l.length&&l.forEach(o=>{s.mergeDateTimeFormat(o,n.datetimeFormats[o])})}if(ks(n.numberFormats)){const l=Object.keys(n.numberFormats);l.length&&l.forEach(o=>{s.mergeNumberFormat(o,n.numberFormats[o])})}}}function tc(s){return ws(Me,null,s,0)}const lc="__INTLIFY_META__",oc=()=>[],ng=()=>!1;let pc=0;function cc(s){return((n,a,e,t)=>s(a,e,ea()||void 0,t))}const ag=()=>{const s=ea();let n=null;return s&&(n=Ri(s)[lc])?{[lc]:n}:null};function Li(s={}){const{__root:n,__injectWithOption:a}=s,e=n===void 0,t=s.flatJson,l=na?os:so;let o=Is(s.inheritLocale)?s.inheritLocale:!0;const p=l(n&&o?n.locale.value:es(s.locale)?s.locale:gt),c=l(n&&o?n.fallbackLocale.value:es(s.fallbackLocale)||qs(s.fallbackLocale)||Cs(s.fallbackLocale)||s.fallbackLocale===!1?s.fallbackLocale:p.value),u=l(Ai(p.value,s)),r=l(Cs(s.datetimeFormats)?s.datetimeFormats:{[p.value]:{}}),i=l(Cs(s.numberFormats)?s.numberFormats:{[p.value]:{}});let d=n?n.missingWarn:Is(s.missingWarn)||ft(s.missingWarn)?s.missingWarn:!0,h=n?n.fallbackWarn:Is(s.fallbackWarn)||ft(s.fallbackWarn)?s.fallbackWarn:!0,y=n?n.fallbackRoot:Is(s.fallbackRoot)?s.fallbackRoot:!0,j=!!s.fallbackFormat,T=Os(s.missing)?s.missing:null,g=Os(s.missing)?cc(s.missing):null,m=Os(s.postTranslation)?s.postTranslation:null,w=n?n.warnHtmlMessage:Is(s.warnHtmlMessage)?s.warnHtmlMessage:!0,f=!!s.escapeParameter;const x=n?n.modifiers:Cs(s.modifiers)?s.modifiers:{};let k=s.pluralRules||n&&n.pluralRules,_;_=(()=>{e&&Kp(null);const L={version:Kf,locale:p.value,fallbackLocale:c.value,messages:u.value,modifiers:x,pluralRules:k,missing:g===null?void 0:g,missingWarn:d,fallbackWarn:h,fallbackFormat:j,unresolving:!0,postTranslation:m===null?void 0:m,warnHtmlMessage:w,escapeParameter:f,messageResolver:s.messageResolver,messageCompiler:s.messageCompiler,__meta:{framework:"vue"}};L.datetimeFormats=r.value,L.numberFormats=i.value,L.__datetimeFormatters=Cs(_)?_.__datetimeFormatters:void 0,L.__numberFormatters=Cs(_)?_.__numberFormatters:void 0,L.__v_emitter=Cs(_)?_.__v_emitter:void 0;const z=Mf(L);return e&&Kp(z),z})(),ne(_,p.value,c.value);function S(){return[p.value,c.value,u.value,r.value,i.value]}const I=ss({get:()=>p.value,set:L=>{_.locale=L,p.value=L}}),J=ss({get:()=>c.value,set:L=>{_.fallbackLocale=L,c.value=L,ne(_,p.value,L)}}),U=ss(()=>u.value),F=ss(()=>Object.keys(u.value).sort()),G=ss(()=>r.value),is=ss(()=>i.value);function ps(){return Os(m)?m:null}function Y(L){m=L,_.postTranslation=L}function us(){return T}function Ns(L){L!==null&&(g=cc(L)),T=L,_.missing=g}function Xs(L,z){return L!=="translate"||!z.resolvedMessage}const Rs=(L,z,ls,vs,Bs,on)=>{S();let Ys;try{e||(_.fallbackContext=n?Lf():void 0),Ys=L(_)}finally{e||(_.fallbackContext=void 0)}if(ls!=="translate exists"&&Ws(Ys)&&Ys===Ot||ls==="translate exists"&&!Ys){const[fn,Fe]=z();if(n&&es(fn)&&Xs(ls,Fe)){y&&(It(h,fn)||Ci(d,fn))&&_a(_o(Va.FALLBACK_TO_ROOT,{key:fn,type:ls}));{const{__v_emitter:an}=_;an&&y&&an.emit("fallback",{type:ls,key:fn,to:"global",groupId:`${ls}:${fn}`})}}return n&&y?vs(n):Bs(fn)}else{if(on(Ys))return Ys;throw Bn(Vs.UNEXPECTED_RETURN_TYPE)}};function zs(...L){return Rs(z=>Reflect.apply(ec,null,[z,...L]),()=>Rl(...L),"translate",z=>Reflect.apply(z.t,z,[...L]),z=>z,z=>es(z))}function ts(...L){const[z,ls,vs]=L;if(vs&&!ks(vs))throw Bn(Vs.INVALID_ARGUMENT);return zs(z,ls,sn({resolvedMessage:!0},vs||{}))}function rs(...L){return Rs(z=>Reflect.apply(Qp,null,[z,...L]),()=>Tl(...L),"datetime format",z=>Reflect.apply(z.d,z,[...L]),()=>bt,z=>es(z)||qs(z))}function fs(...L){return Rs(z=>Reflect.apply(Zp,null,[z,...L]),()=>Al(...L),"number format",z=>Reflect.apply(z.n,z,[...L]),()=>bt,z=>es(z)||qs(z))}function _s(L){return L.map(z=>es(z)||Ws(z)||Is(z)?tc(String(z)):z)}const W={normalize:_s,interpolate:L=>L,type:"vnode"};function H(...L){return Rs(z=>{let ls;const vs=z;try{vs.processor=W,ls=Reflect.apply(ec,null,[vs,...L])}finally{vs.processor=null}return ls},()=>Rl(...L),"translate",z=>z[Ll](...L),z=>[tc(z)],z=>qs(z))}function Q(...L){return Rs(z=>Reflect.apply(Zp,null,[z,...L]),()=>Al(...L),"number format",z=>z[Dl](...L),oc,z=>es(z)||qs(z))}function gs(...L){return Rs(z=>Reflect.apply(Qp,null,[z,...L]),()=>Tl(...L),"datetime format",z=>z[Ml](...L),oc,z=>es(z)||qs(z))}function C(L){k=L,_.pluralRules=k}function P(L,z){return Rs(()=>{if(!L)return!1;const ls=es(z)?z:p.value,vs=V(ls),Bs=_.messageResolver(vs,L);return En(Bs)||xn(Bs)||es(Bs)},()=>[L],"translate exists",ls=>Reflect.apply(ls.te,ls,[L,z]),ng,ls=>Is(ls))}function R(L){let z=null;const ls=gi(_,c.value,p.value);for(let vs=0;vs<ls.length;vs++){const Bs=u.value[ls[vs]]||{},on=_.messageResolver(Bs,L);if(on!=null){z=on;break}}return z}function $(L){const z=R(L);return z??(n?n.tm(L)||{}:{})}function V(L){return u.value[L]||{}}function q(L,z){if(t){const ls={[L]:z};for(const vs in ls)An(ls,vs)&&Ee(ls[vs]);z=ls[L]}u.value[L]=z,_.messages=u.value}function b(L,z){u.value[L]=u.value[L]||{};const ls={[L]:z};if(t)for(const vs in ls)An(ls,vs)&&Ee(ls[vs]);z=ls[L],et(z,u.value[L]),_.messages=u.value}function v(L){return r.value[L]||{}}function A(L,z){r.value[L]=z,_.datetimeFormats=r.value,Jp(_,L,z)}function M(L,z){r.value[L]=sn(r.value[L]||{},z),_.datetimeFormats=r.value,Jp(_,L,z)}function Z(L){return i.value[L]||{}}function K(L,z){i.value[L]=z,_.numberFormats=i.value,sc(_,L,z)}function ns(L,z){i.value[L]=sn(i.value[L]||{},z),_.numberFormats=i.value,sc(_,L,z)}pc++,n&&na&&(Hs(n.locale,L=>{o&&(p.value=L,_.locale=L,ne(_,p.value,c.value))}),Hs(n.fallbackLocale,L=>{o&&(c.value=L,_.fallbackLocale=L,ne(_,p.value,c.value))}));const as={id:pc,locale:I,fallbackLocale:J,get inheritLocale(){return o},set inheritLocale(L){o=L,L&&n&&(p.value=n.locale.value,c.value=n.fallbackLocale.value,ne(_,p.value,c.value))},availableLocales:F,messages:U,get modifiers(){return x},get pluralRules(){return k||{}},get isGlobal(){return e},get missingWarn(){return d},set missingWarn(L){d=L,_.missingWarn=d},get fallbackWarn(){return h},set fallbackWarn(L){h=L,_.fallbackWarn=h},get fallbackRoot(){return y},set fallbackRoot(L){y=L},get fallbackFormat(){return j},set fallbackFormat(L){j=L,_.fallbackFormat=j},get warnHtmlMessage(){return w},set warnHtmlMessage(L){w=L,_.warnHtmlMessage=L},get escapeParameter(){return f},set escapeParameter(L){f=L,_.escapeParameter=L},t:zs,getLocaleMessage:V,setLocaleMessage:q,mergeLocaleMessage:b,getPostTranslationHandler:ps,setPostTranslationHandler:Y,getMissingHandler:us,setMissingHandler:Ns,[Qf]:C};return as.datetimeFormats=G,as.numberFormats=is,as.rt=ts,as.te=P,as.tm=$,as.d=rs,as.n=fs,as.getDateTimeFormat=v,as.setDateTimeFormat=A,as.mergeDateTimeFormat=M,as.getNumberFormat=Z,as.setNumberFormat=K,as.mergeNumberFormat=ns,as[Jf]=a,as[Ll]=H,as[Ml]=gs,as[Dl]=Q,as[Pe]=L=>{_.__v_emitter=L},as[Ol]=()=>{_.__v_emitter=void 0},as}var je=typeof global<"u"?global:typeof self<"u"?self:typeof window<"u"?window:{};function eg(){return Mi().__VUE_DEVTOOLS_GLOBAL_HOOK__}function Mi(){return typeof navigator<"u"&&typeof window<"u"?window:typeof je<"u"?je:{}}const tg=typeof Proxy=="function",lg="devtools-plugin:setup",og="plugin:settings:set";let Oa,Nl;function pg(){var s;return Oa!==void 0||(typeof window<"u"&&window.performance?(Oa=!0,Nl=window.performance):typeof je<"u"&&(!((s=je.perf_hooks)===null||s===void 0)&&s.performance)?(Oa=!0,Nl=je.perf_hooks.performance):Oa=!1),Oa}function cg(){return pg()?Nl.now():Date.now()}class rg{constructor(n,a){this.target=null,this.targetQueue=[],this.onQueue=[],this.plugin=n,this.hook=a;const e={};if(n.settings)for(const o in n.settings){const p=n.settings[o];e[o]=p.defaultValue}const t=`__vue-devtools-plugin-settings__${n.id}`;let l=Object.assign({},e);try{const o=localStorage.getItem(t),p=JSON.parse(o);Object.assign(l,p)}catch{}this.fallbacks={getSettings(){return l},setSettings(o){try{localStorage.setItem(t,JSON.stringify(o))}catch{}l=o},now(){return cg()}},a&&a.on(og,(o,p)=>{o===this.plugin.id&&this.fallbacks.setSettings(p)}),this.proxiedOn=new Proxy({},{get:(o,p)=>this.target?this.target.on[p]:(...c)=>{this.onQueue.push({method:p,args:c})}}),this.proxiedTarget=new Proxy({},{get:(o,p)=>this.target?this.target[p]:p==="on"?this.proxiedOn:Object.keys(this.fallbacks).includes(p)?(...c)=>(this.targetQueue.push({method:p,args:c,resolve:()=>{}}),this.fallbacks[p](...c)):(...c)=>new Promise(u=>{this.targetQueue.push({method:p,args:c,resolve:u})})})}async setRealTarget(n){this.target=n;for(const a of this.onQueue)this.target.on[a.method](...a.args);for(const a of this.targetQueue)a.resolve(await this.target[a.method](...a.args))}}function ig(s,n){const a=s,e=Mi(),t=eg(),l=tg&&a.enableEarlyProxy;if(t&&(e.__VUE_DEVTOOLS_PLUGIN_API_AVAILABLE__||!l))t.emit(lg,s,n);else{const o=l?new rg(a,t):null;(e.__VUE_DEVTOOLS_PLUGINS__=e.__VUE_DEVTOOLS_PLUGINS__||[]).push({pluginDescriptor:a,setupFn:n,proxy:o}),o&&n(o.proxiedTarget)}}const Di="vue-i18n: composer properties",ll={"vue-devtools-plugin-vue-i18n":"Vue I18n DevTools","vue-i18n-resource-inspector":"Vue I18n DevTools","vue-i18n-timeline":"Vue I18n"},ug={"vue-i18n-resource-inspector":"Search for scopes ..."},hg={"vue-i18n-timeline":16764185};let Fl;async function dg(s,n){return new Promise((a,e)=>{try{ig({id:"vue-devtools-plugin-vue-i18n",label:ll["vue-devtools-plugin-vue-i18n"],packageName:"vue-i18n",homepage:"https://vue-i18n.intlify.dev",logo:"https://vue-i18n.intlify.dev/vue-i18n-devtools-logo.png",componentStateTypes:[Di],app:s},t=>{Fl=t,t.on.visitComponentTree(({componentInstance:o,treeNode:p})=>{mg(o,p,n)}),t.on.inspectComponent(({componentInstance:o,instanceData:p})=>{o.vnode.el&&o.vnode.el.__VUE_I18N__&&p&&jg(p,o.vnode.el.__VUE_I18N__)}),t.addInspector({id:"vue-i18n-resource-inspector",label:ll["vue-i18n-resource-inspector"],icon:"language",treeFilterPlaceholder:ug["vue-i18n-resource-inspector"]}),t.on.getInspectorTree(o=>{o.app===s&&o.inspectorId==="vue-i18n-resource-inspector"&&yg(o,n)});const l=new Map;t.on.getInspectorState(async o=>{if(o.app===s&&o.inspectorId==="vue-i18n-resource-inspector")if(t.unhighlightElement(),wg(o,n),o.nodeId==="global"){if(!l.has(o.app)){const[p]=await t.getComponentInstances(o.app);l.set(o.app,p)}t.highlightElement(l.get(o.app))}else{const p=vg(o.nodeId,n);p&&t.highlightElement(p)}}),t.on.editInspectorState(o=>{o.app===s&&o.inspectorId==="vue-i18n-resource-inspector"&&kg(o,n)}),t.addTimelineLayer({id:"vue-i18n-timeline",label:ll["vue-i18n-timeline"],color:hg["vue-i18n-timeline"]}),a(!0)})}catch(t){console.error(t),e(!1)}})}function Oi(s){return s.type.name||s.type.displayName||s.type.__file||"Anonymous"}function mg(s,n,a){const e=a.global;if(s&&s.vnode.el&&s.vnode.el.__VUE_I18N__&&s.vnode.el.__VUE_I18N__!==e){const t={label:`i18n (${Oi(s)} Scope)`,textColor:0,backgroundColor:16764185};n.tags.push(t)}}function jg(s,n){const a=Di;s.state.push({type:a,key:"locale",editable:!0,value:n.locale.value}),s.state.push({type:a,key:"availableLocales",editable:!1,value:n.availableLocales}),s.state.push({type:a,key:"fallbackLocale",editable:!0,value:n.fallbackLocale.value}),s.state.push({type:a,key:"inheritLocale",editable:!0,value:n.inheritLocale}),s.state.push({type:a,key:"messages",editable:!1,value:yo(n.messages.value)}),s.state.push({type:a,key:"datetimeFormats",editable:!1,value:n.datetimeFormats.value}),s.state.push({type:a,key:"numberFormats",editable:!1,value:n.numberFormats.value})}function yo(s){const n={};return Object.keys(s).forEach(a=>{const e=s[a];Os(e)&&"source"in e?n[a]=_g(e):En(e)&&e.loc&&e.loc.source?n[a]=e.loc.source:ks(e)?n[a]=yo(e):n[a]=e}),n}const fg={"<":"&lt;",">":"&gt;",'"':"&quot;","&":"&amp;"};function gg(s){return s.replace(/[<>"&]/g,bg)}function bg(s){return fg[s]||s}function _g(s){return{_custom:{type:"function",display:`<span>ƒ</span> ${s.source?`("${gg(s.source)}")`:"(?)"}`}}}function yg(s,n){s.rootNodes.push({id:"global",label:"Global Scope"});const a=n.global;for(const[e,t]of n.__instances){const l=t;a!==l&&s.rootNodes.push({id:l.id.toString(),label:`${Oi(e)} Scope`})}}function vg(s,n){let a=null;if(s!=="global"){for(const[e,t]of n.__instances.entries())if(t.id.toString()===s){a=e;break}}return a}function Ii(s,n){if(s==="global")return n.global;{const a=Array.from(n.__instances.values()).find(e=>e.id.toString()===s);return a||null}}function wg(s,n){const a=Ii(s.nodeId,n);return a&&(s.state=Cg(a)),null}function Cg(s){const n={},a="Locale related info",e=[{type:a,key:"locale",editable:!0,value:s.locale.value},{type:a,key:"fallbackLocale",editable:!0,value:s.fallbackLocale.value},{type:a,key:"availableLocales",editable:!1,value:s.availableLocales},{type:a,key:"inheritLocale",editable:!0,value:s.inheritLocale}];n[a]=e;const t="Locale messages info",l=[{type:t,key:"messages",editable:!1,value:yo(s.messages.value)}];n[t]=l;{const o="Datetime formats info",p=[{type:o,key:"datetimeFormats",editable:!1,value:s.datetimeFormats.value}];n[o]=p;const c="Datetime formats info",u=[{type:c,key:"numberFormats",editable:!1,value:s.numberFormats.value}];n[c]=u}return n}function $l(s,n){if(Fl){let a;n&&"groupId"in n&&(a=n.groupId,delete n.groupId),Fl.addTimelineEvent({layerId:"vue-i18n-timeline",event:{title:s,groupId:a,time:Date.now(),meta:{},data:n||{},logType:s==="compile-error"?"error":s==="fallback"||s==="missing"?"warning":"default"}})}}function kg(s,n){const a=Ii(s.nodeId,n);if(a){const[e]=s.path;e==="locale"&&es(s.state.value)?a.locale.value=s.state.value:e==="fallbackLocale"&&(es(s.state.value)||qs(s.state.value)||ks(s.state.value))?a.fallbackLocale.value=s.state.value:e==="inheritLocale"&&Is(s.state.value)&&(a.inheritLocale=s.state.value)}}const vo={tag:{type:[String,Object]},locale:{type:String},scope:{type:String,validator:s=>s==="parent"||s==="global",default:"parent"},i18n:{type:Object}};function xg({slots:s},n){return n.length===1&&n[0]==="default"?(s.default?s.default():[]).reduce((e,t)=>[...e,...t.type===bs?t.children:[t]],[]):n.reduce((a,e)=>{const t=s[e];return t&&(a[e]=t()),a},Ts())}function Ni(){return bs}function Pg(s){return qs(s)&&!es(s[0])}function Fi(s,n,a,e){const{slots:t,attrs:l}=n;return()=>{const o={part:!0};let p=Ts();s.locale&&(o.locale=s.locale),es(s.format)?o.key=s.format:ks(s.format)&&(es(s.format.key)&&(o.key=s.format.key),p=Object.keys(s.format).reduce((d,h)=>a.includes(h)?sn(Ts(),d,{[h]:s.format[h]}):d,Ts()));const c=e(s.value,o,p);let u=[o.key];qs(c)?u=c.map((d,h)=>{const y=t[d.type],j=y?y({[d.type]:d.value,index:h,parts:c}):[d.value];return Pg(j)&&(j[0].key=`${d.type}-${h}`),j}):es(c)&&(u=[c]);const r=sn(Ts(),l),i=es(s.tag)||ks(s.tag)?s.tag:Ni();return Oe(i,r,u)}}const Eg=Gs({name:"i18n-d",props:sn({value:{type:[Number,Date],required:!0},format:{type:[String,Object]}},vo),setup(s,n){const a=s.i18n||nn({useScope:s.scope,__useComponent:!0});return Fi(s,n,Pi,(...e)=>a[Ml](...e))}}),rc=Eg,Sg=Gs({name:"i18n-n",props:sn({value:{type:Number,required:!0},format:{type:[String,Object]}},vo),setup(s,n){const a=s.i18n||nn({useScope:s.scope,__useComponent:!0});return Fi(s,n,Ei,(...e)=>a[Dl](...e))}}),ic=Sg,Tg=Gs({name:"i18n-t",props:sn({},{keypath:{type:String,required:!0},plural:{type:[Number,String],validator:s=>Ws(s)||!isNaN(s)}},vo),setup(s,n){const{slots:a,attrs:e}=n,t=s.i18n||nn({useScope:s.scope,__useComponent:!0});return()=>{const l=Object.keys(a).filter(i=>i[0]!=="_"),o=Ts();s.locale&&(o.locale=s.locale),s.plural!==void 0&&(o.plural=es(s.plural)?+s.plural:s.plural);const p=xg(n,l),c=t[Ll](s.keypath,p,o),u=sn(Ts(),e),r=es(s.tag)||ks(s.tag)?s.tag:Ni();return Oe(r,u,c)}}}),uc=Tg;function Ag(s,...n){const a=Cs(n[0])?n[0]:{};(Is(a.globalInstall)?a.globalInstall:!0)&&([uc.name,"I18nT"].forEach(t=>s.component(t,uc)),[ic.name,"I18nN"].forEach(t=>s.component(t,ic)),[rc.name,"I18nD"].forEach(t=>s.component(t,rc)))}const Rg=zn("global-vue-i18n");function Lg(s={}){const n=Is(s.globalInjection)?s.globalInjection:!0,a=new Map,[e,t]=Mg(s),l=zn("vue-i18n");function o(r){return a.get(r)||null}function p(r,i){a.set(r,i)}function c(r){a.delete(r)}const u={async install(r,...i){if(r.__VUE_I18N__=u,r.__VUE_I18N_SYMBOL__=l,r.provide(r.__VUE_I18N_SYMBOL__,u),Cs(i[0])){const y=i[0];u.__composerExtend=y.__composerExtend}let d=null;n&&(d=qg(r,u.global)),Ag(r,...i);const h=r.unmount;r.unmount=()=>{d&&d(),u.dispose(),h()};{if(!await dg(r,u))throw Bn(Vs.CANNOT_SETUP_VUE_DEVTOOLS_PLUGIN);const j=pi(),T=t;T[Pe]&&T[Pe](j),j.on("*",$l)}},get global(){return t},dispose(){e.stop()},__instances:a,__getInstance:o,__setInstance:p,__deleteInstance:c};return u}function nn(s={}){const n=ea();if(n==null)throw Bn(Vs.MUST_BE_CALL_SETUP_TOP);if(!n.isCE&&n.appContext.app!=null&&!n.appContext.app.__VUE_I18N_SYMBOL__)throw Bn(Vs.NOT_INSTALLED);const a=Dg(n),e=Ig(a),t=Ri(n),l=Og(s,t);if(l==="global")return sg(e,s,t),e;if(l==="parent"){let c=Ng(a,n,s.__useComponent);return c==null&&(_a(_o(Va.NOT_FOUND_PARENT_SCOPE)),c=e),c}const o=a;let p=o.__getInstance(n);if(p==null){const c=sn({},s);"__i18n"in t&&(c.__i18n=t.__i18n),e&&(c.__root=e),p=Li(c),o.__composerExtend&&(p[Il]=o.__composerExtend(p)),$g(o,n,p),o.__setInstance(n,p)}else if(l==="local")throw Bn(Vs.DUPLICATE_USE_I18N_CALLING);return p}function Mg(s){const n=Hl(),a=n.run(()=>Li(s));if(a==null)throw Bn(Vs.UNEXPECTED_ERROR);return[n,a]}function Dg(s){const n=jn(s.isCE?Rg:s.appContext.app.__VUE_I18N_SYMBOL__);if(!n)throw Bn(s.isCE?Vs.NOT_INSTALLED_WITH_PROVIDE:Vs.UNEXPECTED_ERROR);return n}function Og(s,n){return Dt(s)?"__i18n"in n?"local":"global":s.useScope?s.useScope:"local"}function Ig(s){return s.global}function Ng(s,n,a=!1){let e=null;const t=n.root;let l=Fg(n,a);for(;l!=null&&(e=s.__getInstance(l),!(e!=null||t===l));)l=l.parent;return e}function Fg(s,n=!1){return s==null?null:n&&s.vnode.ctx||s.parent}function $g(s,n,a){let e=null;wn(()=>{if(n.vnode.el){n.vnode.el.__VUE_I18N__=a,e=pi();const t=a;t[Pe]&&t[Pe](e),e.on("*",$l)}},n),eo(()=>{const t=a;n.vnode.el&&n.vnode.el.__VUE_I18N__&&(e&&e.off("*",$l),t[Ol]&&t[Ol](),delete n.vnode.el.__VUE_I18N__),s.__deleteInstance(n);const l=t[Il];l&&(l(),delete t[Il])},n)}const Bg=["locale","fallbackLocale","availableLocales"],hc=["t","rt","d","n","tm","te"];function qg(s,n){const a=Object.create(null);return Bg.forEach(t=>{const l=Object.getOwnPropertyDescriptor(n,t);if(!l)throw Bn(Vs.UNEXPECTED_ERROR);const o=Fs(l.value)?{get(){return l.value.value},set(p){l.value.value=p}}:{get(){return l.get&&l.get()}};Object.defineProperty(a,t,o)}),s.config.globalProperties.$i18n=a,hc.forEach(t=>{const l=Object.getOwnPropertyDescriptor(n,t);if(!l||!l.value)throw Bn(Vs.UNEXPECTED_ERROR);Object.defineProperty(s.config.globalProperties,`$${t}`,l)}),()=>{delete s.config.globalProperties.$i18n,hc.forEach(t=>{delete s.config.globalProperties[`$${t}`]})}}Ef(tf);Sf(wf);Tf(gi);{const s=wj();s.__INTLIFY__=!0,lf(s.__INTLIFY_DEVTOOLS_GLOBAL_HOOK__)}Xf();const dc=s=>{if(typeof window>"u"||typeof document>"u")return;const n=document.documentElement;if(s==="auto"){const a=window.matchMedia("(prefers-color-scheme: dark)").matches;n.setAttribute("data-bs-theme",a?"dark":"light"),localStorage.removeItem("theme")}else n.setAttribute("data-bs-theme",s),localStorage.setItem("theme",s)},zg={tags:"标签",articles:"文章",words:"字数",prevPage:"上一页",nextPage:"下一页",pagination:"分页导航",search:"搜索",theme:"主题",language:"语言",menu:"菜单",close:"关闭",searchPlaceholder:"搜索文章标题或正文...",searchNoResults:"未找到匹配的文章",searchUnavailable:"搜索索引加载失败",searchEnterHint:"Enter 打开",searchEscHint:"Esc 关闭",searchResultsLabel:"条结果",home:"首页",categories:"分类",resources:"资源",about:"关于",greeting:"你好，我是",greetingPrefix:"//",avatar:"头像",developer:"开发者",wordUnit:"字",recentPosts:"最近文章",notes:"笔记",projects:"项目",topics:"课题",seeMore:"查看更多",countPosts:"{count} 篇",countProjects:"{count} 个项目",countTopics:"{count} 个课题",countWords:"{count} 字",copyCode:"复制代码",copyFailed:"复制失败",copyArticle:"复制文章",copyTable:"复制表格",copied:"已复制",anchorHeading:"置顶当前标题",uncategorized:"未分类",updatedAt:"更新于",postReadingTime:"{minutes} 分钟",articleReadingTime:"阅读约 {minutes} 分钟",tableOfContents:"本页目录",openToc:"打开目录",backToTop:"回到顶部",resourceSubtitle:"生物信息学与结构生物学领域常用工具",resourceItems:"项",experience:"经历",introduction:"自我介绍",contact:"联系我",thanks:"感谢您的关注！",designedByPrefix:"由",designedBySuffix:"设计",clearFilter:"清除筛选",backToArticle:"返回文章"},Vg={tags:"Tags",articles:"Articles",words:"Words",prevPage:"Previous",nextPage:"Next",pagination:"Pagination",search:"Search",theme:"Theme",language:"Language",menu:"Menu",close:"Close",searchPlaceholder:"Search articles by title or content...",searchNoResults:"No matching articles found",searchUnavailable:"Search index failed to load",searchEnterHint:"Enter to open",searchEscHint:"Esc to close",searchResultsLabel:"results",home:"Home",categories:"Categories",resources:"Resources",about:"About",greeting:"Hello, I'm",greetingPrefix:"//",avatar:"Avatar",developer:"Developer",wordUnit:"words",recentPosts:"Recent Posts",notes:"Notes",projects:"Projects",topics:"Topics",seeMore:"See More",countPosts:"{count} posts",countProjects:"{count} projects",countTopics:"{count} topics",countWords:"{count} words",copyCode:"Copy code",copyFailed:"Copy failed",copyArticle:"Copy article",copyTable:"Copy table",copied:"Copied",anchorHeading:"Anchor to heading",uncategorized:"Uncategorized",updatedAt:"Updated at",postReadingTime:"{minutes} min",articleReadingTime:"Reading about {minutes} minutes",tableOfContents:"Table of Contents",openToc:"Open table of contents",backToTop:"Back to top",resourceSubtitle:"Commonly used tools in bioinformatics and structural biology",resourceItems:"items",experience:"Experience",introduction:"Introduction",contact:"Contact Me",thanks:"Thank you for your attention!",designedByPrefix:"Designed by",designedBySuffix:"",clearFilter:"Clear filter",backToArticle:"Back to article"},un={url:"https://zorrooz.github.io",author:"zorrooz",email:"zorrooz@163.com",github:"https://github.com/zorrooz",startYear:2025},Sn=["zh-CN","en-US"],ua=["/zh","/en"],ol=["auto","light","dark"],$i=s=>s==="en"?Sn[1]:Sn[0],Bi=typeof window<"u"&&localStorage.getItem("locale")||Sn[0];typeof document<"u"&&(document.documentElement.lang=Bi);const _t=Lg({locale:Bi,fallbackLocale:Sn[0],messages:{[Sn[0]]:zg,[Sn[1]]:Vg},globalInjection:!0}),wo=fj("app",()=>{const s=os(typeof window<"u"&&localStorage.getItem("theme")||"auto"),n=os(typeof window<"u"&&localStorage.getItem("locale")||Sn[0]);return{theme:s,locale:n,setLocale:o=>{n.value=o,_t.global.locale.value=o,typeof window<"u"&&localStorage.setItem("locale",o),typeof document<"u"&&(document.documentElement.lang=o)},toggleTheme:()=>{const p=(ol.indexOf(s.value)+1)%ol.length;s.value=ol[p],dc(s.value)},initTheme:()=>{dc(s.value)},initLocale:()=>{if(typeof window<"u"){const o=window.location.pathname.match(/^\/(zh|en)(\/|$)/);o&&(n.value=$i(o[1]))}_t.global.locale.value=n.value,typeof document<"u"&&(document.documentElement.lang=n.value)}}}),Qn=s=>`${_t.global.locale.value===Sn[1]?ua[1]:ua[0]}${s.startsWith("/")?s:`/${s}`}`,qi=`Hello, I am zorrooz, currently focused on structural analysis and functional studies of biological macromolecules, while also exploring computational biology.
`,zi=[{year:"2021 – 2025",title:"Lanzhou University · B.S. in Bioinformatics",desc:"Basic programming training and foundational knowledge in biology"},{year:"2025 – present",title:"Lanzhou University · M.S. in Biology",desc:"Structural and functional studies of biological macromolecules"}],Vi=[{title:"Common Tools",items:[{name:"Python · R · Bash · Git",desc:"Primary languages and version control tools for daily research and development"},{name:"Nextflow · Snakemake",desc:"Bioinformatics pipeline orchestration and workflow management"},{name:"RELION · cryoSPARC",desc:"Cryo-EM single-particle reconstruction and 3D structure validation"},{name:"VS Code · Ubuntu · WSL",desc:"Development environment and terminal workflow"}]},{title:"Fields",items:[{name:"Structural Biology",desc:"Cryo-EM"},{name:"Computational Biology",desc:"Binder design and virtual screening"},{name:"Programming",desc:"R · Python · C · Perl · Bash · JavaScript"}]}],Ui=[{label:"Email",value:"zorrooz@163.com",link:"mailto:zorrooz@163.com",icon:"fas fa-envelope"},{label:"GitHub",value:"zorrooz",link:"https://github.com/zorrooz",icon:"fab fa-github"}],Ug={introduction:qi,experience:zi,section:Vi,contacts:Ui},Hg=Object.freeze(Object.defineProperty({__proto__:null,contacts:Ui,default:Ug,experience:zi,introduction:qi,section:Vi},Symbol.toStringTag,{value:"Module"})),Hi=`你好，我是 zorrooz，当前专注于生物大分子的结构解析与功能研究，同时涉猎计算生物学。
`,Gi=[{year:"2021 – 2025",title:"兰州大学 本科 · 生物信息学",desc:"基本编程训练与生物学基础知识"},{year:"2025 – 至今",title:"兰州大学 硕士 · 生物学",desc:"生物大分子的结构与功能研究"}],Wi=[{title:"常用工具",items:[{name:"Python · R · Bash · Git",desc:"日常研究与开发的主力语言与版本控制工具"},{name:"Nextflow · Snakemake",desc:"生信流水线编排与工作流管理"},{name:"RELION · cryoSPARC",desc:"冷冻电镜单颗粒重构与三维结构验证"},{name:"VS Code · Ubuntu · WSL",desc:"开发环境与终端工作流"}]},{title:"领域",items:[{name:"结构生物学",desc:"冷冻电镜"},{name:"计算生物学",desc:"binder 设计与虚拟筛选"},{name:"编程",desc:"R · Python · C · Perl · Bash · JavaScript"}]}],Ki=[{label:"邮箱",value:"zorrooz@163.com",link:"mailto:zorrooz@163.com",icon:"fas fa-envelope"},{label:"GitHub",value:"zorrooz",link:"https://github.com/zorrooz",icon:"fab fa-github"}],Gg={introduction:Hi,experience:Gi,section:Wi,contacts:Ki},Wg=Object.freeze(Object.defineProperty({__proto__:null,contacts:Ki,default:Gg,experience:Gi,introduction:Hi,section:Wi},Symbol.toStringTag,{value:"Module"})),Kg=JSON.parse('[{"title":"Notes","items":[{"name":"ComputationalBiology","title":"Computational Biology","desc":"Mainstream tools for protein design and virtual screening","tags":["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology","Virtual Screening","Molecular Docking","AutoDock Vina"],"stats":{"postsCount":2,"totalWords":1211,"latestDate":"2026-08-04"},"categories":[{"key":"protein-design","title":"Protein Design","articles":[{"title":"Mainstream Tools for Protein Design","articleUrl":"/article/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en","wordCount":618,"tags":["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology"]}],"stats":{"postsCount":1,"totalWords":618,"latestDate":"2026-08-04"}},{"key":"virtual-screening","title":"Virtual Screening","articles":[{"title":"Mainstream Tools for Virtual Screening","articleUrl":"/article/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools-en","wordCount":593,"tags":["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"]}],"stats":{"postsCount":1,"totalWords":593,"latestDate":"2026-08-04"}}],"root":"/article/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en"},{"name":"Programming","title":"Programming","desc":"Detailed tutorials on R, Python, Linux, and Bash programming","tags":["Bash","Shell","Scripting","Tutorial","Linux","Command Line","Python","Advanced","Getting Started","NumPy","Pandas","Data Processing","R Language","ggplot2","Data Visualization","tidyverse","dplyr"],"stats":{"postsCount":8,"totalWords":1972,"latestDate":"2026-08-04"},"categories":[{"key":"python","title":"Python","articles":[{"title":"Python Advanced: Functions, Classes, and Modules","articleUrl":"/article/notes/Programming/python/python-advanced/python-advanced-en","wordCount":244,"tags":["Python","Advanced","Tutorial"]},{"title":"Introduction to Python Programming: Environment, Syntax, and Data Types","articleUrl":"/article/notes/Programming/python/python-basics/python-basics-en","wordCount":333,"tags":["Python","Getting Started","Tutorial"]},{"title":"Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas","articleUrl":"/article/notes/Programming/python/python-data/python-data-en","wordCount":314,"tags":["Python","NumPy","Pandas","Data Processing","Tutorial"]}],"stats":{"postsCount":3,"totalWords":891,"latestDate":"2026-08-04"}},{"key":"r","title":"R","articles":[{"title":"R Language Primer: Data Structures, Vectorization, and Functions","articleUrl":"/article/notes/Programming/r/r-basics/r-basics-en","wordCount":256,"tags":["R Language","Getting Started","Tutorial"]},{"title":"R ggplot2 Data Visualization: Grammar of Graphics and Common Charts","articleUrl":"/article/notes/Programming/r/r-ggplot2/r-ggplot2-en","wordCount":157,"tags":["R Language","ggplot2","Data Visualization","Tutorial"]},{"title":"R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning","articleUrl":"/article/notes/Programming/r/r-tidyverse/r-tidyverse-en","wordCount":245,"tags":["R Language","tidyverse","dplyr","Tutorial"]}],"stats":{"postsCount":3,"totalWords":658,"latestDate":"2026-08-04"}},{"key":"linux","title":"Linux","articles":[{"title":"Linux Command Line Basics: File System, Permissions, and Text Processing","articleUrl":"/article/notes/Programming/linux/linux-basics/linux-basics-en","wordCount":238,"tags":["Linux","Command Line","Tutorial"]}],"stats":{"postsCount":1,"totalWords":238,"latestDate":"2026-08-04"}},{"key":"bash","title":"Bash","articles":[{"title":"Bash Programming: Variables, Conditionals, Loops, and Practical Scripts","articleUrl":"/article/notes/Programming/bash/bash-scripting/bash-scripting-en","wordCount":185,"tags":["Bash","Shell","Scripting","Tutorial"]}],"stats":{"postsCount":1,"totalWords":185,"latestDate":"2026-08-04"}}],"root":"/article/notes/Programming/python/python-advanced/python-advanced-en"},{"name":"StructuralBiology","title":"Structural Biology","desc":"Cryo-EM data processing, visualization of biomacromolecular structures, and atomic modeling","tags":["Cryo-Electron Microscopy","cryo-EM","Review","Data Processing","RELION","Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial","Structure Visualization","PyMOL","ChimeraX"],"stats":{"postsCount":4,"totalWords":2670,"latestDate":"2026-08-04"},"categories":[{"key":"cryoem","title":"Cryo-EM","articles":[{"title":"Review of Cryo-EM Technology","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en","wordCount":977,"tags":["Cryo-Electron Microscopy","cryo-EM","Review"]},{"title":"Cryo-EM Single Particle Analysis: Full Data Processing Workflow","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow-en","wordCount":765,"tags":["Cryo-Electron Microscopy","cryo-EM","Data Processing","RELION"]}],"stats":{"postsCount":2,"totalWords":1742,"latestDate":"2026-08-04"}},{"key":"visualization","title":"Structural Visualization","articles":[{"title":"Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice","articleUrl":"/article/notes/StructuralBiology/visualization/structure-visualization/structure-visualization-en","wordCount":362,"tags":["Structure Visualization","PyMOL","ChimeraX","Tutorial"]}],"stats":{"postsCount":1,"totalWords":362,"latestDate":"2026-08-04"}},{"key":"modeling","title":"Atomic Modeling","articles":[{"title":"Atomic Modeling and Refinement","articleUrl":"/article/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling-en","wordCount":566,"tags":["Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial"]}],"stats":{"postsCount":1,"totalWords":566,"latestDate":"2026-08-04"}}],"root":"/article/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en"}]},{"title":"Projects","items":[{"name":"Plotit","title":"plotit","desc":"Declarative plotting package for R — a verb-prefix API based on ggplot2, in early development stage","tags":["R","ggplot2","Data Visualization","Declarative API"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-12-30"},"categories":[],"root":"","github":"https://github.com/zorrooz/plotit","status":"active","language":"R","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/plotit"},{"name":"ZorroozBlog","title":"zorrooz.github.io","desc":"Personal blog system — built with Vue 3 + Vite, supporting bilingual Chinese/English, theme switching, and Markdown article management","tags":["Vue","Vite","Blog","Markdown"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-08-29"},"categories":[],"root":"","github":"https://github.com/zorrooz/zorrooz.github.io","status":"active","language":"Vue","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/"}]},{"title":"Topics","items":[{"name":"StructureOfProteinDemo","title":"Protein Structure Determination Demo Project","desc":"Protein structure placeholder demo","tags":["Protein","Structure","Cryo-EM"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2026-04-12"},"categories":[],"root":"","doi":"10.1234/demo-structure.2026","status":"completed","journal":"Journal of Demo Science","year":2026,"authors":["Demo Author A","Demo Author B"],"url":"https://www.demo-structure.org"}]}]'),Xg=Object.freeze(Object.defineProperty({__proto__:null,default:Kg},Symbol.toStringTag,{value:"Module"})),Yg=JSON.parse('[{"title":"笔记","items":[{"name":"ComputationalBiology","title":"计算生物学","desc":"蛋白质设计与虚拟筛选的主流工具","tags":["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学","虚拟筛选","分子对接","AutoDock Vina"],"stats":{"postsCount":2,"totalWords":1899,"latestDate":"2026-08-04"},"categories":[{"key":"protein-design","title":"蛋白质设计","articles":[{"title":"蛋白质设计主流工具","articleUrl":"/article/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools","wordCount":1002,"tags":["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学"]}],"stats":{"postsCount":1,"totalWords":1002,"latestDate":"2026-08-04"}},{"key":"virtual-screening","title":"虚拟筛选","articles":[{"title":"虚拟筛选主流工具","articleUrl":"/article/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools","wordCount":897,"tags":["虚拟筛选","分子对接","AutoDock Vina","计算生物学"]}],"stats":{"postsCount":1,"totalWords":897,"latestDate":"2026-08-04"}}],"root":"/article/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools"},{"name":"Programming","title":"编程","desc":"R、Python、Linux 与 Bash 编程的详细教程","tags":["Bash","Shell","脚本编程","教程","Linux","命令行","Python","进阶","入门","NumPy","Pandas","数据处理","R语言","ggplot2","数据可视化","tidyverse","dplyr"],"stats":{"postsCount":8,"totalWords":2903,"latestDate":"2026-08-04"},"categories":[{"key":"python","title":"Python","articles":[{"title":"Python 进阶：函数、类与模块","articleUrl":"/article/notes/Programming/python/python-advanced/python-advanced","wordCount":376,"tags":["Python","进阶","教程"]},{"title":"Python 编程入门：环境、语法与数据类型","articleUrl":"/article/notes/Programming/python/python-basics/python-basics","wordCount":495,"tags":["Python","入门","教程"]},{"title":"Python 数据处理实战：文件 IO、正则与 NumPy/Pandas","articleUrl":"/article/notes/Programming/python/python-data/python-data","wordCount":473,"tags":["Python","NumPy","Pandas","数据处理","教程"]}],"stats":{"postsCount":3,"totalWords":1344,"latestDate":"2026-08-04"}},{"key":"r","title":"R","articles":[{"title":"R 语言入门：数据结构、向量化与函数","articleUrl":"/article/notes/Programming/r/r-basics/r-basics","wordCount":403,"tags":["R语言","入门","教程"]},{"title":"R ggplot2 数据可视化：图层语法与常用图表","articleUrl":"/article/notes/Programming/r/r-ggplot2/r-ggplot2","wordCount":245,"tags":["R语言","ggplot2","数据可视化","教程"]},{"title":"R tidyverse 数据操作：dplyr 管道与数据清洗","articleUrl":"/article/notes/Programming/r/r-tidyverse/r-tidyverse","wordCount":324,"tags":["R语言","tidyverse","dplyr","教程"]}],"stats":{"postsCount":3,"totalWords":972,"latestDate":"2026-08-04"}},{"key":"linux","title":"Linux","articles":[{"title":"Linux 命令行基础：文件系统、权限与文本处理","articleUrl":"/article/notes/Programming/linux/linux-basics/linux-basics","wordCount":328,"tags":["Linux","命令行","教程"]}],"stats":{"postsCount":1,"totalWords":328,"latestDate":"2026-08-04"}},{"key":"bash","title":"Bash","articles":[{"title":"Bash 编程：变量、条件、循环与实用脚本","articleUrl":"/article/notes/Programming/bash/bash-scripting/bash-scripting","wordCount":259,"tags":["Bash","Shell","脚本编程","教程"]}],"stats":{"postsCount":1,"totalWords":259,"latestDate":"2026-08-04"}}],"root":"/article/notes/Programming/python/python-advanced/python-advanced"},{"name":"StructuralBiology","title":"结构生物学","desc":"冷冻电镜数据处理、生物大分子结构可视化与原子建模","tags":["冷冻电镜","cryo-EM","综述","数据处理","RELION","原子建模","Coot","phenix","结构精修","教程","结构可视化","PyMOL","ChimeraX"],"stats":{"postsCount":4,"totalWords":3960,"latestDate":"2026-08-04"},"categories":[{"key":"cryoem","title":"冷冻电镜","articles":[{"title":"冷冻电镜技术综述","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview","wordCount":1485,"tags":["冷冻电镜","cryo-EM","综述"]},{"title":"冷冻电镜单颗粒分析：数据处理全流程","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow","wordCount":1111,"tags":["冷冻电镜","cryo-EM","数据处理","RELION"]}],"stats":{"postsCount":2,"totalWords":2596,"latestDate":"2026-08-04"}},{"key":"visualization","title":"结构可视化","articles":[{"title":"生物大分子结构可视化：PyMOL 与 ChimeraX 实战","articleUrl":"/article/notes/StructuralBiology/visualization/structure-visualization/structure-visualization","wordCount":529,"tags":["结构可视化","PyMOL","ChimeraX","教程"]}],"stats":{"postsCount":1,"totalWords":529,"latestDate":"2026-08-04"}},{"key":"modeling","title":"原子建模","articles":[{"title":"原子建模与精修","articleUrl":"/article/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling","wordCount":835,"tags":["原子建模","Coot","phenix","结构精修","教程"]}],"stats":{"postsCount":1,"totalWords":835,"latestDate":"2026-08-04"}}],"root":"/article/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview"}]},{"title":"项目","items":[{"name":"Plotit","title":"plotit","desc":"R 声明式绘图包——基于 ggplot2 的动词前缀 API，早期开发阶段","tags":["R","ggplot2","数据可视化","声明式API"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-12-30"},"categories":[],"root":"","github":"https://github.com/zorrooz/plotit","status":"active","language":"R","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/plotit"},{"name":"ZorroozBlog","title":"zorrooz.github.io","desc":"个人博客系统——基于 Vue 3 + Vite 构建，支持中英双语、主题切换与 Markdown 文章管理","tags":["Vue","Vite","博客","Markdown"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-08-29"},"categories":[],"root":"","github":"https://github.com/zorrooz/zorrooz.github.io","status":"active","language":"Vue","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/"}]},{"title":"课题","items":[{"name":"StructureOfProteinDemo","title":"蛋白质结构解析演示课题","desc":"蛋白质结构占位demo","tags":["Protein","Structure","Cryo-EM"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2026-04-12"},"categories":[],"root":"","doi":"10.1234/demo-structure.2026","status":"completed","journal":"Journal of Demo Science","year":2026,"authors":["Demo Author A","Demo Author B"],"url":"https://www.demo-structure.org"}]}]'),Qg=Object.freeze(Object.defineProperty({__proto__:null,default:Yg},Symbol.toStringTag,{value:"Module"})),Jg=[{title:"Mainstream Tools for Protein Design",date:"2026-08-04",tags:["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology"],draft:!1,description:"A panorama of protein design tools: Rosetta energy optimization, AlphaFold structure prediction, RFdiffusion generative design, ESMFold fast prediction, and binder design workflow",relativePath:"ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en",wordCount:618,tagCount:5},{title:"Mainstream Tools for Virtual Screening",date:"2026-08-04",tags:["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"],draft:!1,description:"A comprehensive overview of virtual screening tools: docking software such as AutoDock Vina, Glide, and DOCK, receptor/ligand preparation, scoring functions, high-throughput screening workflows, and validation strategies",relativePath:"ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools-en",wordCount:593,tagCount:4},{title:"Bash Programming: Variables, Conditionals, Loops, and Practical Scripts",date:"2026-08-04",tags:["Bash","Shell","Scripting","Tutorial"],draft:!1,description:"A complete Bash scripting tutorial: variables and parameters, conditional statements, loops, functions, arrays, and string processing, with practical scripts for batch processing and bioinformatics",relativePath:"Programming/bash/bash-scripting/bash-scripting-en",wordCount:185,tagCount:4},{title:"Linux Command Line Basics: File System, Permissions, and Text Processing",date:"2026-08-04",tags:["Linux","Command Line","Tutorial"],draft:!1,description:"Core skills for the Linux command line: file system navigation, file operations, permission management, the three musketeers of text processing (grep/sed/awk), and process management",relativePath:"Programming/linux/linux-basics/linux-basics-en",wordCount:238,tagCount:3},{title:"Python Advanced: Functions, Classes, and Modules",date:"2026-08-04",tags:["Python","Advanced","Tutorial"],draft:!1,description:"Deep dive into Python functional and object-oriented programming: parameter passing, lambda, decorators, classes and inheritance, exception handling, modules and packages",relativePath:"Programming/python/python-advanced/python-advanced-en",wordCount:244,tagCount:3},{title:"Introduction to Python Programming: Environment, Syntax, and Data Types",date:"2026-08-04",tags:["Python","Getting Started","Tutorial"],draft:!1,description:"Learn Python from scratch: environment setup, basic syntax, built-in data types, and control flow, with runnable example code.",relativePath:"Programming/python/python-basics/python-basics-en",wordCount:333,tagCount:3},{title:"Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas",date:"2026-08-04",tags:["Python","NumPy","Pandas","Data Processing","Tutorial"],draft:!1,description:"A complete practical tutorial on file I/O, regular expressions, and the scientific computing stack (NumPy/Pandas), covering common data scenarios in bioinformatics",relativePath:"Programming/python/python-data/python-data-en",wordCount:314,tagCount:5},{title:"R Language Primer: Data Structures, Vectorization, and Functions",date:"2026-08-04",tags:["R Language","Getting Started","Tutorial"],draft:!1,description:"Core R fundamentals: environment setup, vectors/matrices/data frames/lists, vectorized operations, control flow, and function definitions",relativePath:"Programming/r/r-basics/r-basics-en",wordCount:256,tagCount:3},{title:"R ggplot2 Data Visualization: Grammar of Graphics and Common Charts",date:"2026-08-04",tags:["R Language","ggplot2","Data Visualization","Tutorial"],draft:!1,description:"Introduction to the ggplot2 Grammar of Graphics: scatter plots, boxplots, histograms, bar charts, and theme customization for publication-ready figures",relativePath:"Programming/r/r-ggplot2/r-ggplot2-en",wordCount:157,tagCount:4},{title:"R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning",date:"2026-08-04",tags:["R Language","tidyverse","dplyr","Tutorial"],draft:!1,description:"Chain dplyr verbs (filter/select/mutate/group_by/summarise/join) with pipes, and use tidyr to complete a full data cleaning workflow.",relativePath:"Programming/r/r-tidyverse/r-tidyverse-en",wordCount:245,tagCount:4},{title:"Review of Cryo-EM Technology",date:"2026-08-04",tags:["Cryo-Electron Microscopy","cryo-EM","Review"],draft:!1,description:"The technological revolution of single-particle cryo-EM: direct electron detectors, the stability revolution, and the advent of the AI era, along with the future landscape of structural biology",relativePath:"StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en",wordCount:977,tagCount:3},{title:"Cryo-EM Single Particle Analysis: Full Data Processing Workflow",date:"2026-08-04",tags:["Cryo-Electron Microscopy","cryo-EM","Data Processing","RELION"],draft:!1,description:"From micrographs to atomic models: the complete workflow of cryo-EM single particle analysis (SPA), covering preprocessing, particle picking, 2D/3D classification, refinement, and resolution assessment",relativePath:"StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow-en",wordCount:765,tagCount:4},{title:"Atomic Modeling and Refinement",date:"2026-08-04",tags:["Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial"],draft:!1,description:"Complete workflow for atomic model building and refinement from cryo-EM density maps: automatic modeling, Coot manual correction, phenix refinement, model quality assessment metrics",relativePath:"StructuralBiology/modeling/atomic-modeling/atomic-modeling-en",wordCount:566,tagCount:5},{title:"Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice",date:"2026-08-04",tags:["Structure Visualization","PyMOL","ChimeraX","Tutorial"],draft:!1,description:"Visualizing protein/nucleic acid structures with PyMOL and UCSF ChimeraX: PDB data retrieval, rendering modes, script batch processing, and cryo-EM density map display",relativePath:"StructuralBiology/visualization/structure-visualization/structure-visualization-en",wordCount:362,tagCount:4}],Zg=Object.freeze(Object.defineProperty({__proto__:null,default:Jg},Symbol.toStringTag,{value:"Module"})),sb=[{title:"蛋白质设计主流工具",date:"2026-08-04",tags:["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学"],draft:!1,description:"蛋白质设计工具全景：Rosetta 能量优化、AlphaFold 结构预测、RFdiffusion 生成式设计、ESMFold 快速预测与 binder 设计工作流",relativePath:"ComputationalBiology/protein-design/protein-design-tools/protein-design-tools",wordCount:1002,tagCount:5},{title:"虚拟筛选主流工具",date:"2026-08-04",tags:["虚拟筛选","分子对接","AutoDock Vina","计算生物学"],draft:!1,description:"虚拟筛选工具全景：AutoDock Vina、Glide、DOCK 等对接软件，受体/配体准备、打分函数、高通量筛选流程与验证策略",relativePath:"ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools",wordCount:897,tagCount:4},{title:"Bash 编程：变量、条件、循环与实用脚本",date:"2026-08-04",tags:["Bash","Shell","脚本编程","教程"],draft:!1,description:"Bash 脚本编程完整教程：变量与参数、条件判断、循环、函数、数组与字符串处理，附批处理与生信实用脚本",relativePath:"Programming/bash/bash-scripting/bash-scripting",wordCount:259,tagCount:4},{title:"Linux 命令行基础：文件系统、权限与文本处理",date:"2026-08-04",tags:["Linux","命令行","教程"],draft:!1,description:"Linux 命令行核心技能：文件系统导航、文件操作、权限管理、文本处理三剑客（grep/sed/awk）与进程管理",relativePath:"Programming/linux/linux-basics/linux-basics",wordCount:328,tagCount:3},{title:"Python 进阶：函数、类与模块",date:"2026-08-04",tags:["Python","进阶","教程"],draft:!1,description:"深入 Python 函数式与面向对象编程：参数传递、lambda、装饰器、类与继承、异常处理、模块与包",relativePath:"Programming/python/python-advanced/python-advanced",wordCount:376,tagCount:3},{title:"Python 编程入门：环境、语法与数据类型",date:"2026-08-04",tags:["Python","入门","教程"],draft:!1,description:"从零开始学习 Python：环境搭建、基础语法、内置数据类型与流程控制，附带可运行的示例代码",relativePath:"Programming/python/python-basics/python-basics",wordCount:495,tagCount:3},{title:"Python 数据处理实战：文件 IO、正则与 NumPy/Pandas",date:"2026-08-04",tags:["Python","NumPy","Pandas","数据处理","教程"],draft:!1,description:"文件读写、正则表达式与科学计算三件套（NumPy/Pandas）的完整实战教程，覆盖生物信息学常见数据场景",relativePath:"Programming/python/python-data/python-data",wordCount:473,tagCount:5},{title:"R 语言入门：数据结构、向量化与函数",date:"2026-08-04",tags:["R语言","入门","教程"],draft:!1,description:"R 语言核心基础：环境配置、向量/矩阵/数据框/列表、向量化运算、控制流与函数定义",relativePath:"Programming/r/r-basics/r-basics",wordCount:403,tagCount:3},{title:"R ggplot2 数据可视化：图层语法与常用图表",date:"2026-08-04",tags:["R语言","ggplot2","数据可视化","教程"],draft:!1,description:"ggplot2 图层语法入门：散点图、箱线图、直方图、柱状图与主题美化，输出出版级图表",relativePath:"Programming/r/r-ggplot2/r-ggplot2",wordCount:245,tagCount:4},{title:"R tidyverse 数据操作：dplyr 管道与数据清洗",date:"2026-08-04",tags:["R语言","tidyverse","dplyr","教程"],draft:!1,description:"用管道符串联 dplyr 动词（filter/select/mutate/group_by/summarise/join），配合 tidyr 完成完整的数据清洗流程",relativePath:"Programming/r/r-tidyverse/r-tidyverse",wordCount:324,tagCount:4},{title:"冷冻电镜技术综述",date:"2026-08-04",tags:["冷冻电镜","cryo-EM","综述"],draft:!1,description:"冷冻电镜单颗粒分析的技术革命：直接电子探测器、稳定性革命与 AI 时代的到来，以及结构生物学的未来图景",relativePath:"StructuralBiology/cryoem/cryoem-overview/cryoem-overview",wordCount:1485,tagCount:3},{title:"冷冻电镜单颗粒分析：数据处理全流程",date:"2026-08-04",tags:["冷冻电镜","cryo-EM","数据处理","RELION"],draft:!1,description:"从显微照片到原子模型：冷冻电镜单颗粒分析（SPA）完整流程，涵盖预处理、颗粒挑选、2D/3D 分类、精修与分辨率评估",relativePath:"StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow",wordCount:1111,tagCount:4},{title:"原子建模与精修",date:"2026-08-04",tags:["原子建模","Coot","phenix","结构精修","教程"],draft:!1,description:"冷冻电镜密度图的原子模型搭建与精修全流程：自动建模、Coot 手动修正、phenix 精修、模型质量评估指标",relativePath:"StructuralBiology/modeling/atomic-modeling/atomic-modeling",wordCount:835,tagCount:5},{title:"生物大分子结构可视化：PyMOL 与 ChimeraX 实战",date:"2026-08-04",tags:["结构可视化","PyMOL","ChimeraX","教程"],draft:!1,description:"用 PyMOL 与 UCSF ChimeraX 可视化蛋白质/核酸结构：PDB 数据获取、渲染模式、脚本批处理与冷冻电镜密度图展示",relativePath:"StructuralBiology/visualization/structure-visualization/structure-visualization",wordCount:529,tagCount:4}],nb=Object.freeze(Object.defineProperty({__proto__:null,default:sb},Symbol.toStringTag,{value:"Module"})),ab=[{id:1,no:1,title:"Mainstream Tools for Protein Design",date:"2026-08-04",category:["ComputationalBiology","protein-design"],tags:["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology"],preview:"A panorama of protein design tools: Rosetta energy optimization, AlphaFold structure prediction, RFdiffusion generative design, ESMFold fast prediction, and binder design workflow",wordCount:618},{id:2,no:2,title:"Mainstream Tools for Virtual Screening",date:"2026-08-04",category:["ComputationalBiology","virtual-screening"],tags:["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"],preview:"A comprehensive overview of virtual screening tools: docking software such as AutoDock Vina, Glide, and DOCK, receptor/ligand preparation, scoring functions, high-throughput screening workflows, and validation strategies",wordCount:593},{id:3,no:3,title:"Bash Programming: Variables, Conditionals, Loops, and Practical Scripts",date:"2026-08-04",category:["Programming","bash"],tags:["Bash","Shell","Scripting","Tutorial"],preview:"A complete Bash scripting tutorial: variables and parameters, conditional statements, loops, functions, arrays, and string processing, with practical scripts for batch processing and bioinformatics",wordCount:185},{id:4,no:4,title:"Linux Command Line Basics: File System, Permissions, and Text Processing",date:"2026-08-04",category:["Programming","linux"],tags:["Linux","Command Line","Tutorial"],preview:"Core skills for the Linux command line: file system navigation, file operations, permission management, the three musketeers of text processing (grep/sed/awk), and process management",wordCount:238},{id:5,no:5,title:"Python Advanced: Functions, Classes, and Modules",date:"2026-08-04",category:["Programming","python"],tags:["Python","Advanced","Tutorial"],preview:"Deep dive into Python functional and object-oriented programming: parameter passing, lambda, decorators, classes and inheritance, exception handling, modules and packages",wordCount:244},{id:6,no:6,title:"Introduction to Python Programming: Environment, Syntax, and Data Types",date:"2026-08-04",category:["Programming","python"],tags:["Python","Getting Started","Tutorial"],preview:"Learn Python from scratch: environment setup, basic syntax, built-in data types, and control flow, with runnable example code.",wordCount:333},{id:7,no:7,title:"Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas",date:"2026-08-04",category:["Programming","python"],tags:["Python","NumPy","Pandas","Data Processing","Tutorial"],preview:"A complete practical tutorial on file I/O, regular expressions, and the scientific computing stack (NumPy/Pandas), covering common data scenarios in bioinformatics",wordCount:314},{id:8,no:8,title:"R Language Primer: Data Structures, Vectorization, and Functions",date:"2026-08-04",category:["Programming","r"],tags:["R Language","Getting Started","Tutorial"],preview:"Core R fundamentals: environment setup, vectors/matrices/data frames/lists, vectorized operations, control flow, and function definitions",wordCount:256},{id:9,no:9,title:"R ggplot2 Data Visualization: Grammar of Graphics and Common Charts",date:"2026-08-04",category:["Programming","r"],tags:["R Language","ggplot2","Data Visualization","Tutorial"],preview:"Introduction to the ggplot2 Grammar of Graphics: scatter plots, boxplots, histograms, bar charts, and theme customization for publication-ready figures",wordCount:157},{id:10,no:10,title:"R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning",date:"2026-08-04",category:["Programming","r"],tags:["R Language","tidyverse","dplyr","Tutorial"],preview:"Chain dplyr verbs (filter/select/mutate/group_by/summarise/join) with pipes, and use tidyr to complete a full data cleaning workflow.",wordCount:245},{id:11,no:11,title:"Review of Cryo-EM Technology",date:"2026-08-04",category:["StructuralBiology","cryoem"],tags:["Cryo-Electron Microscopy","cryo-EM","Review"],preview:"The technological revolution of single-particle cryo-EM: direct electron detectors, the stability revolution, and the advent of the AI era, along with the future landscape of structural biology",wordCount:977},{id:12,no:12,title:"Cryo-EM Single Particle Analysis: Full Data Processing Workflow",date:"2026-08-04",category:["StructuralBiology","cryoem"],tags:["Cryo-Electron Microscopy","cryo-EM","Data Processing","RELION"],preview:"From micrographs to atomic models: the complete workflow of cryo-EM single particle analysis (SPA), covering preprocessing, particle picking, 2D/3D classification, refinement, and resolution assessment",wordCount:765},{id:13,no:13,title:"Atomic Modeling and Refinement",date:"2026-08-04",category:["StructuralBiology","modeling"],tags:["Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial"],preview:"Complete workflow for atomic model building and refinement from cryo-EM density maps: automatic modeling, Coot manual correction, phenix refinement, model quality assessment metrics",wordCount:566},{id:14,no:14,title:"Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice",date:"2026-08-04",category:["StructuralBiology","visualization"],tags:["Structure Visualization","PyMOL","ChimeraX","Tutorial"],preview:"Visualizing protein/nucleic acid structures with PyMOL and UCSF ChimeraX: PDB data retrieval, rendering modes, script batch processing, and cryo-EM density map display",wordCount:362}],eb=Object.freeze(Object.defineProperty({__proto__:null,default:ab},Symbol.toStringTag,{value:"Module"})),tb=[{id:1,no:1,title:"蛋白质设计主流工具",date:"2026-08-04",category:["计算生物学","protein-design"],tags:["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学"],preview:"蛋白质设计工具全景：Rosetta 能量优化、AlphaFold 结构预测、RFdiffusion 生成式设计、ESMFold 快速预测与 binder 设计工作流",wordCount:1002},{id:2,no:2,title:"虚拟筛选主流工具",date:"2026-08-04",category:["计算生物学","virtual-screening"],tags:["虚拟筛选","分子对接","AutoDock Vina","计算生物学"],preview:"虚拟筛选工具全景：AutoDock Vina、Glide、DOCK 等对接软件，受体/配体准备、打分函数、高通量筛选流程与验证策略",wordCount:897},{id:3,no:3,title:"Bash 编程：变量、条件、循环与实用脚本",date:"2026-08-04",category:["编程","bash"],tags:["Bash","Shell","脚本编程","教程"],preview:"Bash 脚本编程完整教程：变量与参数、条件判断、循环、函数、数组与字符串处理，附批处理与生信实用脚本",wordCount:259},{id:4,no:4,title:"Linux 命令行基础：文件系统、权限与文本处理",date:"2026-08-04",category:["编程","linux"],tags:["Linux","命令行","教程"],preview:"Linux 命令行核心技能：文件系统导航、文件操作、权限管理、文本处理三剑客（grep/sed/awk）与进程管理",wordCount:328},{id:5,no:5,title:"Python 进阶：函数、类与模块",date:"2026-08-04",category:["编程","python"],tags:["Python","进阶","教程"],preview:"深入 Python 函数式与面向对象编程：参数传递、lambda、装饰器、类与继承、异常处理、模块与包",wordCount:376},{id:6,no:6,title:"Python 编程入门：环境、语法与数据类型",date:"2026-08-04",category:["编程","python"],tags:["Python","入门","教程"],preview:"从零开始学习 Python：环境搭建、基础语法、内置数据类型与流程控制，附带可运行的示例代码",wordCount:495},{id:7,no:7,title:"Python 数据处理实战：文件 IO、正则与 NumPy/Pandas",date:"2026-08-04",category:["编程","python"],tags:["Python","NumPy","Pandas","数据处理","教程"],preview:"文件读写、正则表达式与科学计算三件套（NumPy/Pandas）的完整实战教程，覆盖生物信息学常见数据场景",wordCount:473},{id:8,no:8,title:"R 语言入门：数据结构、向量化与函数",date:"2026-08-04",category:["编程","r"],tags:["R语言","入门","教程"],preview:"R 语言核心基础：环境配置、向量/矩阵/数据框/列表、向量化运算、控制流与函数定义",wordCount:403},{id:9,no:9,title:"R ggplot2 数据可视化：图层语法与常用图表",date:"2026-08-04",category:["编程","r"],tags:["R语言","ggplot2","数据可视化","教程"],preview:"ggplot2 图层语法入门：散点图、箱线图、直方图、柱状图与主题美化，输出出版级图表",wordCount:245},{id:10,no:10,title:"R tidyverse 数据操作：dplyr 管道与数据清洗",date:"2026-08-04",category:["编程","r"],tags:["R语言","tidyverse","dplyr","教程"],preview:"用管道符串联 dplyr 动词（filter/select/mutate/group_by/summarise/join），配合 tidyr 完成完整的数据清洗流程",wordCount:324},{id:11,no:11,title:"冷冻电镜技术综述",date:"2026-08-04",category:["结构生物学","cryoem"],tags:["冷冻电镜","cryo-EM","综述"],preview:"冷冻电镜单颗粒分析的技术革命：直接电子探测器、稳定性革命与 AI 时代的到来，以及结构生物学的未来图景",wordCount:1485},{id:12,no:12,title:"冷冻电镜单颗粒分析：数据处理全流程",date:"2026-08-04",category:["结构生物学","cryoem"],tags:["冷冻电镜","cryo-EM","数据处理","RELION"],preview:"从显微照片到原子模型：冷冻电镜单颗粒分析（SPA）完整流程，涵盖预处理、颗粒挑选、2D/3D 分类、精修与分辨率评估",wordCount:1111},{id:13,no:13,title:"原子建模与精修",date:"2026-08-04",category:["结构生物学","modeling"],tags:["原子建模","Coot","phenix","结构精修","教程"],preview:"冷冻电镜密度图的原子模型搭建与精修全流程：自动建模、Coot 手动修正、phenix 精修、模型质量评估指标",wordCount:835},{id:14,no:14,title:"生物大分子结构可视化：PyMOL 与 ChimeraX 实战",date:"2026-08-04",category:["结构生物学","visualization"],tags:["结构可视化","PyMOL","ChimeraX","教程"],preview:"用 PyMOL 与 UCSF ChimeraX 可视化蛋白质/核酸结构：PDB 数据获取、渲染模式、脚本批处理与冷冻电镜密度图展示",wordCount:529}],lb=Object.freeze(Object.defineProperty({__proto__:null,default:tb},Symbol.toStringTag,{value:"Module"})),ob=[{name:"Plotit",title:"plotit",desc:"Declarative plotting package for R — a verb-prefix API based on ggplot2, in early development stage",date:"2025-12-30",tags:["R","ggplot2","Data Visualization","Declarative API"],github:"https://github.com/zorrooz/plotit",url:"https://zorrooz.github.io/plotit",status:"active",language:"R",stars:0,license:"MIT",version:""},{name:"ZorroozBlog",title:"zorrooz.github.io",desc:"Personal blog system — built with Vue 3 + Vite, supporting bilingual Chinese/English, theme switching, and Markdown article management",date:"2025-08-29",tags:["Vue","Vite","Blog","Markdown"],github:"https://github.com/zorrooz/zorrooz.github.io",url:"https://zorrooz.github.io/",status:"active",language:"Vue",stars:0,license:"MIT",version:""}],pb=Object.freeze(Object.defineProperty({__proto__:null,default:ob},Symbol.toStringTag,{value:"Module"})),cb=[{name:"Plotit",title:"plotit",desc:"R 声明式绘图包——基于 ggplot2 的动词前缀 API，早期开发阶段",date:"2025-12-30",tags:["R","ggplot2","数据可视化","声明式API"],github:"https://github.com/zorrooz/plotit",url:"https://zorrooz.github.io/plotit",status:"active",language:"R",stars:0,license:"MIT",version:""},{name:"ZorroozBlog",title:"zorrooz.github.io",desc:"个人博客系统——基于 Vue 3 + Vite 构建，支持中英双语、主题切换与 Markdown 文章管理",date:"2025-08-29",tags:["Vue","Vite","博客","Markdown"],github:"https://github.com/zorrooz/zorrooz.github.io",url:"https://zorrooz.github.io/",status:"active",language:"Vue",stars:0,license:"MIT",version:""}],rb=Object.freeze(Object.defineProperty({__proto__:null,default:cb},Symbol.toStringTag,{value:"Module"})),ib=JSON.parse(`[{"title":"Data Analysis","children":[{"title":"Data Exploration","children":[{"title":"R Language","items":[{"name":"tidyverse","url":"https://www.tidyverse.org/","desc":"A suite of R packages for data science, covering data import, cleaning, transformation, and visualization"},{"name":"data.table","url":"https://rdatatable.gitlab.io/data.table/","desc":"High-performance data manipulation package, extremely fast for processing data frames with millions of rows"},{"name":"readr","url":"https://readr.tidyverse.org/","desc":"A tool for fast reading of CSV/TSV and other tabular files, with automatic column type inference"}]},{"title":"Python","items":[{"name":"NumPy","url":"https://numpy.org/","desc":"The fundamental library for scientific computing in Python, providing multidimensional arrays and vectorized operations"},{"name":"Pandas","url":"https://pandas.pydata.org/","desc":"A powerful data analysis library based on Python, providing efficient data manipulation and processing capabilities"},{"name":"Polars","url":"https://www.pola.rs/","desc":"A high-performance data processing library based on Apache Arrow, supporting multiple language interfaces"},{"name":"Jupyter","url":"https://jupyter.org/","desc":"An interactive notebook environment that supports mixing code, text, and visualizations"}]}]},{"title":"Data Visualization","children":[{"title":"R Language","items":[{"name":"ggplot2","url":"https://ggplot2.tidyverse.org/","desc":"A powerful data visualization package in R, based on the grammar of graphics"},{"name":"plotly (R)","url":"https://plotly.com/r/","desc":"The R interface to the interactive chart library, capable of generating web-interactive graphics"},{"name":"patchwork","url":"https://patchwork.data-imaginist.com/","desc":"A multi-panel composition tool that easily combines ggplot graphics using + and / symbols"},{"name":"RColorBrewer","url":"https://cran.r-project.org/web/packages/RColorBrewer/","desc":"Publication-quality color schemes, providing carefully designed discrete palettes"}]},{"title":"Python","items":[{"name":"Matplotlib","url":"https://matplotlib.org/","desc":"A widely used plotting library in Python, supporting many chart types"},{"name":"Seaborn","url":"https://seaborn.pydata.org/","desc":"A statistical visualization library based on Matplotlib, with built-in attractive themes and statistical plots"},{"name":"Plotly","url":"https://plotly.com/python/","desc":"An interactive visualization library supporting scatter, 3D, maps, and many other chart types"},{"name":"Altair","url":"https://altair-viz.github.io/","desc":"A declarative statistical visualization library based on Vega-Lite"}]}]},{"title":"Statistical Analysis","children":[{"title":"R Language","items":[{"name":"R Stats","url":"https://www.r-project.org/","desc":"The built-in statistical analysis functionality in R, covering a wide range of statistical methods and models"},{"name":"rstatix","url":"https://rpkgs.datanovia.com/rstatix/","desc":"Pipe-friendly wrappers for statistical tests (t-tests, ANOVA, rank-sum tests, etc.)"},{"name":"broom","url":"https://broom.tidymodels.org/","desc":"A tool for tidying statistical model outputs into clean data frames"},{"name":"lme4","url":"https://cran.r-project.org/web/packages/lme4/","desc":"Linear and nonlinear mixed-effects models, suitable for repeated-measures data"}]},{"title":"Python","items":[{"name":"SciPy Stats","url":"https://scipy.org/","desc":"The statistics module in the Python scientific computing library SciPy, providing many statistical distributions and test methods"},{"name":"statsmodels","url":"https://www.statsmodels.org/","desc":"A statistical modeling library supporting regression, time series, hypothesis testing, and more"},{"name":"pingouin","url":"https://pingouin-stats.org/","desc":"A simple and user-friendly statistical testing library covering common parametric and nonparametric tests"}]}]},{"title":"Machine Learning","children":[{"title":"Python","items":[{"name":"scikit-learn","url":"https://scikit-learn.org/","desc":"A machine learning library based on Python, providing a rich set of algorithms and tools"},{"name":"PyTorch","url":"https://pytorch.org/","desc":"A deep learning framework supporting dynamic computation graphs and efficient tensor operations"},{"name":"TensorFlow","url":"https://www.tensorflow.org/","desc":"A deep learning framework with a mature ecosystem, supporting production-level deployment"},{"name":"XGBoost","url":"https://xgboost.readthedocs.io/","desc":"A gradient boosting tree library, the top choice for tabular data competitions and engineering"},{"name":"LightGBM","url":"https://lightgbm.readthedocs.io/","desc":"A gradient boosting framework from Microsoft, with fast training speed and low memory usage"}]},{"title":"R Language","items":[{"name":"tidymodels","url":"https://www.tidymodels.org/","desc":"A unified modeling framework in R, covering data preprocessing, modeling, and evaluation"}]}]}]},{"title":"Omics","children":[{"title":"Data Platforms","children":[{"title":"Sequence Databases","items":[{"name":"NCBI","url":"https://www.ncbi.nlm.nih.gov/","desc":"The National Center for Biotechnology Information in the United States, providing important databases such as GenBank"},{"name":"EMBL-EBI","url":"https://www.ebi.ac.uk/","desc":"The European Molecular Biology Laboratory, providing a variety of bioinformatics resources and tools"},{"name":"DDBJ","url":"https://www.ddbj.nig.ac.jp/","desc":"The DNA Data Bank of Japan, providing storage and access for nucleic acid and protein sequence data"},{"name":"UniProt","url":"https://www.uniprot.org/","desc":"The most comprehensive database for protein sequences and functional annotation"},{"name":"Ensembl","url":"https://www.ensembl.org/","desc":"A genome annotation and comparative genomics database for vertebrates"},{"name":"UCSC Genome Browser","url":"https://genome.ucsc.edu/","desc":"A genome browser supporting visualization of multiple species genomes and custom tracks"}]},{"title":"Sequencing Data","items":[{"name":"SRA","url":"https://www.ncbi.nlm.nih.gov/sra","desc":"The NCBI Sequence Read Archive, storing raw high-throughput sequencing data"},{"name":"GEO","url":"https://www.ncbi.nlm.nih.gov/geo/","desc":"A gene expression database containing microarray and sequencing expression data"},{"name":"ENA","url":"https://www.ebi.ac.uk/ena","desc":"The European Nucleotide Archive, Europe's storage center for raw sequencing data"}]},{"title":"Protein Interactions and Pathways","items":[{"name":"STRING","url":"https://string-db.org/","desc":"A protein-protein interaction database supporting multiple species"},{"name":"KEGG","url":"https://www.kegg.jp/","desc":"A comprehensive database of pathways, genes, and compounds"},{"name":"GO","url":"http://geneontology.org/","desc":"Gene Ontology, providing a standardized system for gene function annotation"}]}]},{"title":"Analysis Software","children":[{"title":"Command Line","items":[{"name":"FastQC","url":"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/","desc":"A quality assessment tool for sequencing data, generating visual QC reports"},{"name":"Trimmomatic","url":"http://www.usadellab.org/cms/?page=trimmomatic","desc":"A read trimming tool for sequencing data, removing adapters and low-quality bases"},{"name":"cutadapt","url":"https://cutadapt.readthedocs.io/","desc":"A command-line adapter removal tool supporting multiple sequencing platforms"},{"name":"BWA","url":"https://github.com/lh3/bwa","desc":"A short-read alignment tool, the standard choice for genome alignment"},{"name":"STAR","url":"https://github.com/alexdobin/STAR","desc":"An ultra-fast splice-aware aligner for RNA-seq"},{"name":"HISAT2","url":"https://daehwankimlab.github.io/hisat2/","desc":"A graph-based index alignment tool for RNA-seq"},{"name":"SAMtools","url":"https://www.htslib.org/","desc":"A toolset for processing SAM/BAM files, the standard for post-alignment processing"},{"name":"bcftools","url":"https://samtools.github.io/bcftools/","desc":"A tool for processing VCF/BCF variant files, supporting filtering and statistics"},{"name":"GATK","url":"https://gatk.broadinstitute.org/","desc":"The Genome Analysis Toolkit, the industry standard for variant detection"},{"name":"MultiQC","url":"https://multiqc.info/","desc":"Aggregates multiple QC reports into a single interactive HTML report"}]},{"title":"Python","items":[{"name":"Biopython","url":"https://biopython.org/","desc":"The core Python library for bioinformatics, covering sequences, structures, and file parsing"},{"name":"scikit-bio","url":"https://scikit-bio.org/","desc":"A Python library for bioinformatics, providing various tools for biological data processing and analysis"},{"name":"scanpy","url":"https://scanpy.readthedocs.io/","desc":"A Python library for single-cell transcriptomics analysis, deeply integrated with AnnData"}]},{"title":"R Language","items":[{"name":"Bioconductor","url":"https://bioconductor.org/","desc":"A collection of bioinformatics software packages based on R, covering genomics, transcriptomics, and other fields"},{"name":"Seurat","url":"https://satijalab.org/seurat/","desc":"A mainstream R package for single-cell transcriptomics analysis, providing a complete analysis workflow"},{"name":"DESeq2","url":"https://bioconductor.org/packages/release/bioc/html/DESeq2.html","desc":"The standard tool for differential expression analysis based on the negative binomial model"},{"name":"edgeR","url":"https://bioconductor.org/packages/release/bioc/html/edgeR.html","desc":"A differential expression analysis tool based on empirical Bayes"},{"name":"limma","url":"https://bioconductor.org/packages/release/bioc/html/limma.html","desc":"A differential analysis tool using linear models, applicable to both microarray and sequencing data"}]}]}]},{"title":"Cryo-EM Structure Determination","children":[{"title":"Software Tools","children":[{"title":"3D Reconstruction","items":[{"name":"RELION","url":"https://www3.mrc-lmb.cam.ac.uk/groups/scheres/index.php","desc":"A Bayesian approach-based software for cryo-EM 3D reconstruction, widely used for high-resolution structure determination"},{"name":"cryoSPARC","url":"https://cryoem.slac.stanford.edu/","desc":"A high-performance cryo-EM data processing platform with an intuitive user interface"},{"name":"cisTEM","url":"https://cistem.org/","desc":"A free, open-source cryo-EM processing software supporting a complete single-particle workflow"},{"name":"SPHIRE","url":"https://sphire.mpg.de/","desc":"A modular cryo-EM processing suite supporting high-resolution reconstruction"}]},{"title":"Preprocessing","items":[{"name":"MotionCor2","url":"https://emcore.ucsf.edu/ucsf-software","desc":"A frame-by-frame motion correction tool that eliminates sample drift caused by the electron beam"},{"name":"CTFFIND4","url":"https://grigoriefflab.umassmed.edu/ctffind4","desc":"A classic tool for CTF estimation, evaluating defocus and astigmatism"},{"name":"Gctf","url":"https://www.mrc-lmb.cam.ac.uk/kzhang/Gctf/","desc":"A GPU-accelerated CTF estimation tool with high speed and precision"},{"name":"Warp","url":"http://www.warpem.com/","desc":"A fast preprocessing tool integrating motion correction, CTF estimation, and particle picking"}]},{"title":"Particle Picking","items":[{"name":"crYOLO","url":"https://cryolo.readthedocs.io/","desc":"A deep learning-based particle picking tool supporting real-time picking"},{"name":"Topaz","url":"https://topaz.readthedocs.io/","desc":"A CNN-based particle picking tool with excellent performance on low signal-to-noise ratio samples"}]},{"title":"Modeling and Visualization","items":[{"name":"PyMOL","url":"https://pymol.org/","desc":"A molecular visualization software, the standard tool for scripted rendering and structure analysis"},{"name":"UCSF ChimeraX","url":"https://www.cgl.ucsf.edu/chimerax/","desc":"A modern molecular visualization software with outstanding density map processing capabilities"},{"name":"Coot","url":"https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/","desc":"A manual model building tool for interactive adjustment of density maps and models"},{"name":"phenix","url":"https://phenix-online.org/","desc":"A structure refinement and validation suite supporting cryo-EM and crystallography"},{"name":"ISOLDE","url":"https://isolde.cimr.cam.ac.uk/","desc":"An interactive real-time molecular dynamics modeling tool based on ChimeraX"},{"name":"ModelAngelo","url":"https://github.com/3dem/model-angelo","desc":"A deep learning-based automatic atomic modeling tool for density maps"}]}]},{"title":"Database Resources","children":[{"title":"Structures and Density Maps","items":[{"name":"EMDB","url":"https://www.ebi.ac.uk/emdb/","desc":"The Electron Microscopy Data Bank, storing and distributing 3D density maps from cryo-EM"},{"name":"EMPIAR","url":"https://www.ebi.ac.uk/empiar/","desc":"The Electron Microscopy Public Image Archive, providing raw movies and particle data"},{"name":"RCSB PDB","url":"https://www.rcsb.org/","desc":"The Protein Data Bank, containing biological macromolecular structures determined by X-ray, NMR, and cryo-EM"},{"name":"wwPDB","url":"https://www.wwpdb.org/","desc":"The worldwide PDB coordination body, uniformly managing structural data standards"},{"name":"PDBe","url":"https://www.ebi.ac.uk/pdbe/","desc":"The European PDB node, providing structural data and visualization tools"},{"name":"AlphaFold DB","url":"https://alphafold.ebi.ac.uk/","desc":"The AlphaFold protein structure database, covering hundreds of millions of proteins"}]}]}]}]`),ub=Object.freeze(Object.defineProperty({__proto__:null,default:ib},Symbol.toStringTag,{value:"Module"})),hb=JSON.parse('[{"title":"数据分析","children":[{"title":"数据探索","children":[{"title":"R 语言","items":[{"name":"tidyverse","url":"https://www.tidyverse.org/","desc":"R 语言中用于数据科学的一整套包，涵盖数据导入、清洗、转换和可视化"},{"name":"data.table","url":"https://rdatatable.gitlab.io/data.table/","desc":"高性能数据操作包，百万行级数据框处理速度极快"},{"name":"readr","url":"https://readr.tidyverse.org/","desc":"快速读取 CSV/TSV 等表格文件的工具，自动推断列类型"}]},{"title":"Python","items":[{"name":"NumPy","url":"https://numpy.org/","desc":"Python 科学计算基础库，提供多维数组与向量化运算"},{"name":"Pandas","url":"https://pandas.pydata.org/","desc":"基于 Python 的强大数据分析库，提供高效的数据操作和处理功能"},{"name":"Polars","url":"https://www.pola.rs/","desc":"基于 Apache Arrow 的高性能数据处理库，支持多语言接口"},{"name":"Jupyter","url":"https://jupyter.org/","desc":"交互式笔记本环境，支持代码、文本与可视化混排"}]}]},{"title":"数据可视化","children":[{"title":"R 语言","items":[{"name":"ggplot2","url":"https://ggplot2.tidyverse.org/","desc":"R 语言中功能强大的数据可视化包，基于语法图形学"},{"name":"plotly（R）","url":"https://plotly.com/r/","desc":"交互式图表库的 R 接口，可生成网页端可交互图形"},{"name":"patchwork","url":"https://patchwork.data-imaginist.com/","desc":"多图拼接工具，用 + 和 / 符号轻松组合 ggplot 图形"},{"name":"RColorBrewer","url":"https://cran.r-project.org/web/packages/RColorBrewer/","desc":"出版级配色方案，提供精心设计的离散色板"}]},{"title":"Python","items":[{"name":"Matplotlib","url":"https://matplotlib.org/","desc":"Python 中广泛使用的绘图库，支持多种图表类型"},{"name":"Seaborn","url":"https://seaborn.pydata.org/","desc":"基于 Matplotlib 的统计可视化库，内置美观主题与统计图"},{"name":"Plotly","url":"https://plotly.com/python/","desc":"交互式可视化库，支持散点、3D、地图等多种图表"},{"name":"Altair","url":"https://altair-viz.github.io/","desc":"基于 Vega-Lite 的声明式统计可视化库"}]}]},{"title":"统计分析","children":[{"title":"R 语言","items":[{"name":"R Stats","url":"https://www.r-project.org/","desc":"R 语言内置的统计分析功能，涵盖广泛的统计方法和模型"},{"name":"rstatix","url":"https://rpkgs.datanovia.com/rstatix/","desc":"管道友好的统计检验封装（t 检验、方差分析、秩和检验等）"},{"name":"broom","url":"https://broom.tidymodels.org/","desc":"把统计模型输出整理为整洁数据框的工具"},{"name":"lme4","url":"https://cran.r-project.org/web/packages/lme4/","desc":"线性与非线性混合效应模型，适用于重复测量数据"}]},{"title":"Python","items":[{"name":"SciPy Stats","url":"https://scipy.org/","desc":"Python 科学计算库 SciPy 中的统计模块，提供多种统计分布和检验方法"},{"name":"statsmodels","url":"https://www.statsmodels.org/","desc":"统计建模库，支持回归、时间序列、假设检验等"},{"name":"pingouin","url":"https://pingouin-stats.org/","desc":"简洁易用的统计检验库，覆盖常用参数与非参数检验"}]}]},{"title":"机器学习","children":[{"title":"Python","items":[{"name":"scikit-learn","url":"https://scikit-learn.org/","desc":"基于 Python 的机器学习库，提供丰富的算法和工具"},{"name":"PyTorch","url":"https://pytorch.org/","desc":"深度学习框架，支持动态计算图和高效的张量运算"},{"name":"TensorFlow","url":"https://www.tensorflow.org/","desc":"深度学习框架，生态成熟，支持生产级部署"},{"name":"XGBoost","url":"https://xgboost.readthedocs.io/","desc":"梯度提升树库，表格数据竞赛与工程的首选"},{"name":"LightGBM","url":"https://lightgbm.readthedocs.io/","desc":"微软出品的梯度提升框架，训练速度快、内存占用低"}]},{"title":"R 语言","items":[{"name":"tidymodels","url":"https://www.tidymodels.org/","desc":"R 语言统一的建模框架，覆盖数据预处理、建模与评估"}]}]}]},{"title":"组学","children":[{"title":"数据平台","children":[{"title":"序列数据库","items":[{"name":"NCBI","url":"https://www.ncbi.nlm.nih.gov/","desc":"美国国家生物技术信息中心，提供 GenBank 等重要数据库"},{"name":"EMBL-EBI","url":"https://www.ebi.ac.uk/","desc":"欧洲分子生物学实验室，提供多种生物信息学资源和工具"},{"name":"DDBJ","url":"https://www.ddbj.nig.ac.jp/","desc":"日本 DNA 数据库，提供核酸、蛋白序列数据的存储和访问"},{"name":"UniProt","url":"https://www.uniprot.org/","desc":"最全面的蛋白质序列与功能注释数据库"},{"name":"Ensembl","url":"https://www.ensembl.org/","desc":"脊椎动物基因组注释与比较基因组学数据库"},{"name":"UCSC Genome Browser","url":"https://genome.ucsc.edu/","desc":"基因组浏览器，支持多物种基因组可视化与自定义轨道"}]},{"title":"测序数据","items":[{"name":"SRA","url":"https://www.ncbi.nlm.nih.gov/sra","desc":"NCBI 原始测序数据仓库，存储高通量测序原始数据"},{"name":"GEO","url":"https://www.ncbi.nlm.nih.gov/geo/","desc":"基因表达数据库，收录芯片与测序表达数据"},{"name":"ENA","url":"https://www.ebi.ac.uk/ena","desc":"欧洲核酸档案库，欧洲的原始测序数据存储中心"}]},{"title":"蛋白互作与通路","items":[{"name":"STRING","url":"https://string-db.org/","desc":"蛋白-蛋白相互作用数据库，支持多种物种"},{"name":"KEGG","url":"https://www.kegg.jp/","desc":"通路、基因与化合物综合数据库"},{"name":"GO","url":"http://geneontology.org/","desc":"基因本体论，提供标准化的基因功能注释体系"}]}]},{"title":"分析软件","children":[{"title":"命令行","items":[{"name":"FastQC","url":"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/","desc":"测序数据质量评估工具，生成可视化质控报告"},{"name":"Trimmomatic","url":"http://www.usadellab.org/cms/?page=trimmomatic","desc":"测序读段修剪工具，去除接头与低质量碱基"},{"name":"cutadapt","url":"https://cutadapt.readthedocs.io/","desc":"命令行接头去除工具，支持多种测序平台"},{"name":"BWA","url":"https://github.com/lh3/bwa","desc":"短读段比对工具，基因组比对的标准选择"},{"name":"STAR","url":"https://github.com/alexdobin/STAR","desc":"RNA-seq 剪接感知比对器，速度极快"},{"name":"HISAT2","url":"https://daehwankimlab.github.io/hisat2/","desc":"基于图索引的 RNA-seq 比对工具"},{"name":"SAMtools","url":"https://www.htslib.org/","desc":"SAM/BAM 文件处理工具集，比对后处理的标准"},{"name":"bcftools","url":"https://samtools.github.io/bcftools/","desc":"VCF/BCF 变异文件处理工具，支持过滤与统计"},{"name":"GATK","url":"https://gatk.broadinstitute.org/","desc":"基因组分析工具包，变异检测的行业标准"},{"name":"MultiQC","url":"https://multiqc.info/","desc":"把多个质控报告汇总为一个交互式 HTML 报告"}]},{"title":"Python","items":[{"name":"Biopython","url":"https://biopython.org/","desc":"生物信息学核心 Python 库，涵盖序列、结构与文件解析"},{"name":"scikit-bio","url":"https://scikit-bio.org/","desc":"用于生物信息学的 Python 库，提供多种生物数据处理和分析工具"},{"name":"scanpy","url":"https://scanpy.readthedocs.io/","desc":"单细胞转录组分析 Python 库，与 AnnData 深度集成"}]},{"title":"R 语言","items":[{"name":"Bioconductor","url":"https://bioconductor.org/","desc":"基于 R 语言的生物信息学软件包集合，涵盖基因组学、转录组学等领域"},{"name":"Seurat","url":"https://satijalab.org/seurat/","desc":"单细胞转录组分析主流 R 包，提供完整分析流程"},{"name":"DESeq2","url":"https://bioconductor.org/packages/release/bioc/html/DESeq2.html","desc":"基于负二项模型的差异表达分析标准工具"},{"name":"edgeR","url":"https://bioconductor.org/packages/release/bioc/html/edgeR.html","desc":"基于经验贝叶斯的差异表达分析工具"},{"name":"limma","url":"https://bioconductor.org/packages/release/bioc/html/limma.html","desc":"线性模型差异分析工具，芯片与测序数据通用"}]}]}]},{"title":"电镜结构解析","children":[{"title":"软件工具","children":[{"title":"三维重构","items":[{"name":"RELION","url":"https://www3.mrc-lmb.cam.ac.uk/groups/scheres/index.php","desc":"基于贝叶斯方法的冷冻电镜三维重构软件，广泛应用于高分辨率结构解析"},{"name":"cryoSPARC","url":"https://cryoem.slac.stanford.edu/","desc":"高性能的冷冻电镜数据处理平台，提供直观的用户界面"},{"name":"cisTEM","url":"https://cistem.org/","desc":"免费开源的冷冻电镜处理软件，支持完整的单颗粒流程"},{"name":"SPHIRE","url":"https://sphire.mpg.de/","desc":"模块化冷冻电镜处理套件，支持高分辨率重构"}]},{"title":"预处理","items":[{"name":"MotionCor2","url":"https://emcore.ucsf.edu/ucsf-software","desc":"逐帧运动校正工具，消除电子束导致的样品漂移"},{"name":"CTFFIND4","url":"https://grigoriefflab.umassmed.edu/ctffind4","desc":"CTF 估计经典工具，评估欠焦与像散"},{"name":"Gctf","url":"https://www.mrc-lmb.cam.ac.uk/kzhang/Gctf/","desc":"GPU 加速的 CTF 估计工具，速度快精度高"},{"name":"Warp","url":"http://www.warpem.com/","desc":"集运动校正、CTF 与颗粒挑选于一体的快速预处理工具"}]},{"title":"颗粒挑选","items":[{"name":"crYOLO","url":"https://cryolo.readthedocs.io/","desc":"基于深度学习的颗粒挑选工具，支持实时挑选"},{"name":"Topaz","url":"https://topaz.readthedocs.io/","desc":"基于 CNN 的颗粒挑选工具，对低信噪比样品表现优异"}]},{"title":"建模与可视化","items":[{"name":"PyMOL","url":"https://pymol.org/","desc":"分子可视化软件，脚本化渲染与结构分析的标准工具"},{"name":"UCSF ChimeraX","url":"https://www.cgl.ucsf.edu/chimerax/","desc":"现代分子可视化软件，密度图处理能力突出"},{"name":"Coot","url":"https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/","desc":"手动模型修正工具，密度图与模型交互调整"},{"name":"phenix","url":"https://phenix-online.org/","desc":"结构精修与验证套件，支持冷冻电镜与晶体学"},{"name":"ISOLDE","url":"https://isolde.cimr.cam.ac.uk/","desc":"交互式实时分子动力学建模工具，基于 ChimeraX"},{"name":"ModelAngelo","url":"https://github.com/3dem/model-angelo","desc":"基于深度学习的密度图自动原子建模工具"}]}]},{"title":"数据库资源","children":[{"title":"结构与密度图","items":[{"name":"EMDB","url":"https://www.ebi.ac.uk/emdb/","desc":"电子显微镜数据银行，存储和分发冷冻电镜三维密度图"},{"name":"EMPIAR","url":"https://www.ebi.ac.uk/empiar/","desc":"冷冻电镜原始数据档案，提供原始电影与颗粒数据"},{"name":"RCSB PDB","url":"https://www.rcsb.org/","desc":"蛋白质数据银行，包含通过X射线、NMR和冷冻电镜确定的生物大分子结构"},{"name":"wwPDB","url":"https://www.wwpdb.org/","desc":"全球 PDB 协调机构，统一管理结构数据标准"},{"name":"PDBe","url":"https://www.ebi.ac.uk/pdbe/","desc":"欧洲 PDB 节点，提供结构数据与可视化工具"},{"name":"AlphaFold DB","url":"https://alphafold.ebi.ac.uk/","desc":"AlphaFold 预测蛋白质结构数据库，覆盖数亿蛋白质"}]}]}]}]'),db=Object.freeze(Object.defineProperty({__proto__:null,default:hb},Symbol.toStringTag,{value:"Module"})),mb=[{name:"Advanced",count:1},{name:"AlphaFold",count:1},{name:"Atomic Modeling",count:1},{name:"AutoDock Vina",count:1},{name:"Bash",count:1},{name:"ChimeraX",count:1},{name:"Command Line",count:1},{name:"Computational Biology",count:2},{name:"Coot",count:1},{name:"Cryo-Electron Microscopy",count:2},{name:"cryo-EM",count:2},{name:"Data Processing",count:2},{name:"Data Visualization",count:1},{name:"dplyr",count:1},{name:"Getting Started",count:2},{name:"ggplot2",count:1},{name:"Linux",count:1},{name:"Molecular Docking",count:1},{name:"NumPy",count:1},{name:"Pandas",count:1},{name:"phenix",count:1},{name:"Protein Design",count:1},{name:"PyMOL",count:1},{name:"Python",count:3},{name:"R Language",count:3},{name:"RELION",count:1},{name:"Review",count:1},{name:"RFdiffusion",count:1},{name:"Rosetta",count:1},{name:"Scripting",count:1},{name:"Shell",count:1},{name:"Structure Refinement",count:1},{name:"Structure Visualization",count:1},{name:"tidyverse",count:1},{name:"Tutorial",count:10},{name:"Virtual Screening",count:1}],jb=Object.freeze(Object.defineProperty({__proto__:null,default:mb},Symbol.toStringTag,{value:"Module"})),fb=[{name:"蛋白质设计",count:1},{name:"分子对接",count:1},{name:"计算生物学",count:2},{name:"脚本编程",count:1},{name:"教程",count:10},{name:"结构精修",count:1},{name:"结构可视化",count:1},{name:"进阶",count:1},{name:"冷冻电镜",count:2},{name:"命令行",count:1},{name:"入门",count:2},{name:"数据处理",count:2},{name:"数据可视化",count:1},{name:"虚拟筛选",count:1},{name:"原子建模",count:1},{name:"综述",count:1},{name:"AlphaFold",count:1},{name:"AutoDock Vina",count:1},{name:"Bash",count:1},{name:"ChimeraX",count:1},{name:"Coot",count:1},{name:"cryo-EM",count:2},{name:"dplyr",count:1},{name:"ggplot2",count:1},{name:"Linux",count:1},{name:"NumPy",count:1},{name:"Pandas",count:1},{name:"phenix",count:1},{name:"PyMOL",count:1},{name:"Python",count:3},{name:"R语言",count:3},{name:"RELION",count:1},{name:"RFdiffusion",count:1},{name:"Rosetta",count:1},{name:"Shell",count:1},{name:"tidyverse",count:1}],gb=Object.freeze(Object.defineProperty({__proto__:null,default:fb},Symbol.toStringTag,{value:"Module"})),bb=[{name:"StructureOfProteinDemo",title:"Protein Structure Determination Demo Project",desc:"Protein structure placeholder demo",date:"2026-04-12",tags:["Protein","Structure","Cryo-EM"],doi:"10.1234/demo-structure.2026",url:"https://www.demo-structure.org",status:"completed",journal:"Journal of Demo Science",year:2026,authors:["Demo Author A","Demo Author B"]}],_b=Object.freeze(Object.defineProperty({__proto__:null,default:bb},Symbol.toStringTag,{value:"Module"})),yb=[{name:"StructureOfProteinDemo",title:"蛋白质结构解析演示课题",desc:"蛋白质结构占位demo",date:"2026-04-12",tags:["Protein","Structure","Cryo-EM"],doi:"10.1234/demo-structure.2026",url:"https://www.demo-structure.org",status:"completed",journal:"Journal of Demo Science",year:2026,authors:["Demo Author A","Demo Author B"]}],vb=Object.freeze(Object.defineProperty({__proto__:null,default:yb},Symbol.toStringTag,{value:"Module"})),wb=`<h1>Mainstream Tools for Protein Design</h1>
<p>Protein design is an important direction in computational biology: given a target function, design an amino acid sequence that folds into a specific structure. This article reviews the mainstream tool spectrum from classical energy functions to generative AI.</p>
<h2>1. Two Paradigms of Design Problems</h2>
<h3>1.1 Sequence Design</h3>
<p>Given a target structure (backbone), find the most stable amino acid sequence:</p>
<pre><code>structure (backbone) → sequence design → amino acid sequence → experimental validation
</code></pre>
<h3>1.2 Structure Design (De Novo Design)</h3>
<p>Design new structures/functions from scratch:</p>
<pre><code>functional requirement → structure generation → sequence design → experimental validation
</code></pre>
<h2>2. Rosetta: The Energy Function Paradigm</h2>
<p>Rosetta is the "veteran empire" in structure prediction and design, based on a hybrid energy function of physics and knowledge.</p>
<h3>2.1 Core Energy Function</h3>
<ul>
<li><strong>REF2015</strong>: the classic all-atom energy function (van der Waals, hydrogen bonding, solvation, electrostatics, and dihedral angle preferences)</li>
<li>Score = weighted sum of energy terms; lower is better</li>
<li>Rosetta's "design" is essentially <strong>energy minimization</strong>: searching for low-energy states in sequence/structure space</li>
</ul>
<h3>2.2 Common Protocols</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Sequence design (FixBB: fixed backbone design)</span>
rosetta_scripts.default.linuxgccrelease \\
  -s input.pdb \\
  -parser:protocol fixbb.xml \\
  -out:prefix designed_

<span class="hljs-comment"># Core snippet of fixbb.xml</span>
<span class="hljs-comment"># &#x3C;TaskOperations></span>
<span class="hljs-comment">#   &#x3C;RestrictToRepacking name="repack"/></span>
<span class="hljs-comment">#   &#x3C;OperateOnResidueSubset name="no_design_core" selector="..." /></span>
<span class="hljs-comment"># &#x3C;/TaskOperations></span>
<span class="hljs-comment"># &#x3C;RosettaScripts></span>
<span class="hljs-comment">#   &#x3C;PackRotamersMover name="pack" task_operations="repack"/></span>
<span class="hljs-comment"># &#x3C;/RosettaScripts></span>
</code></pre>
<h3>2.3 Classic Applications of Rosetta</h3>
<ul>
<li><strong>Enzyme design</strong>: design of scaffolds for catalytic active sites (e.g., Kemp eliminase)</li>
<li><strong>Protein-protein interface design</strong>: engineering binding specificity</li>
<li><strong>Binder design</strong>: high-affinity binding proteins</li>
<li><strong>Thermostability engineering</strong>: computational mutations to increase Tm</li>
</ul>
<h3>2.4 Advantages and Limitations</h3>
<table>
<thead>
<tr>
<th>Advantages</th>
<th>Limitations</th>
</tr>
</thead>
<tbody>
<tr>
<td>Physically interpretable, finely controllable</td>
<td>Computationally intensive, requires empirical parameter tuning</td>
</tr>
<tr>
<td>Supports arbitrary backbones and modifications</td>
<td>Limited prediction for flexible regions</td>
</tr>
<tr>
<td>Mature ecosystem (ample documentation, tutorials)</td>
<td>Steep learning curve</td>
</tr>
</tbody>
</table>
<h2>3. AlphaFold: Prediction as a Design Aid</h2>
<h3>3.1 AlphaFold2 (2021)</h3>
<ul>
<li>End-to-end deep learning: sequence → structure (MSA + attention mechanism)</li>
<li>Accuracy: atomic-level accuracy at CASP14 (GDT > 90)</li>
<li>Value for design: <strong>rapidly verify the foldability of designed sequences</strong></li>
</ul>
<h3>3.2 AlphaFold3 (2024)</h3>
<ul>
<li>Unified framework for predicting protein, nucleic acid, small molecule, and ion complexes</li>
<li>Introduces diffusion models and pairwise representations</li>
<li>Significance: designed targets (binder–ligand complex structures) can be directly predicted</li>
</ul>
<pre><code class="hljs language-bash"><span class="hljs-comment"># ColabFold: a lightweight implementation of AlphaFold2 (GPU-friendly)</span>
colabfold_batch input.fasta out_dir \\
  --num-recycle 3 --model-type alphafold2_multimer_v3
</code></pre>
<h3>3.3 Key Usage: Validating Designs</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Confidence indicators for designed sequences</span>
<span class="hljs-comment"># pLDDT > 90: high confidence</span>
<span class="hljs-comment"># PAE &#x3C; 5: inter-domain relative positions trustworthy</span>
<span class="hljs-comment"># Design validation pipeline: sequence → AF2 prediction → comparison with target structure (RMSD)</span>
</code></pre>
<h2>4. Generative Design: RFdiffusion and ProteinMPNN</h2>
<h3>4.1 RFdiffusion (2023, Baker Lab)</h3>
<p>Structure generator based on <strong>denoising diffusion models</strong>:</p>
<ul>
<li>Input: functional motifs, shape constraints, binding targets</li>
<li>Output: novel protein backbones</li>
<li>Applications: binder design, symmetric oligomers, enzyme active-site scaffolds</li>
</ul>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RFdiffusion example: designing a binding protein</span>
run_inference.py \\
  scaffoldguided.scaffoldguided=True \\
  scaff_loader.contigs_map.json=contigs.json \\
  inference.output_prefix=outputs/design \\
  inference.num_designs=100 \\
  potentials.guidance_scale=3.0
</code></pre>
<h3>4.2 ProteinMPNN (2022)</h3>
<p><strong>Sequence design neural network</strong>: given a backbone → sequence, forming a "generate-encode" loop with RFdiffusion:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># ProteinMPNN sequence design</span>
python protein_mpnn_run.py \\
  --pdbpath scaffold.pdb \\
  --out_folder mpnn_out \\
  --num_seq_per_target 8 \\
  --batch_size 4
</code></pre>
<p><strong>Standard workflow (Baker Lab paradigm)</strong>:</p>
<pre><code>① RFdiffusion: generate a backbone
   ↓
② ProteinMPNN: backbone → multiple candidate sequences (sampling temperature 0.1-0.3)
   ↓
③ AF2 reverse validation: sequence → predicted structure → RMSD comparison with backbone
   ↓
④ Screen candidates with high pLDDT / low RMSD
   ↓
⑤ Experimental expression validation (yeast/phage display, ELISA)
</code></pre>
<h3>4.3 Other Generative Tools</h3>
<ul>
<li><strong>Chroma</strong>: another diffusion model (generation + sequence design in one)</li>
<li><strong>Genie</strong>: a faster protein diffusion model</li>
<li><strong>ESMFold</strong> (Meta, 2022): ultra-fast prediction without MSA (&#x3C;1 second per sequence), suitable for large-scale design screening</li>
<li><strong>ESM3</strong>: multimodal generation (sequence + structure + function)</li>
</ul>
<h2>5. Special Topic: Binder Design</h2>
<p>Binder (binding protein) design is currently the most active application scenario:</p>
<h3>5.1 Classic Workflow (Combining Existing Tools)</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. Hotspot analysis on target surface (Rosetta)</span>
<span class="hljs-comment"># 2. Hotspot-guided binder backbone generation (RFdiffusion / or hotspot constraints)</span>
<span class="hljs-comment"># 3. ProteinMPNN sequence design</span>
<span class="hljs-comment"># 4. AF2 complex prediction validation (binder + target docking)</span>
<span class="hljs-comment"># 5. Experimental screening (yeast display, FACS sorting)</span>
</code></pre>
<h3>5.2 Key Considerations</h3>
<ul>
<li>Determination of binding hotspot residues: experimental mutation scanning (alanine scan) or computational (Rosetta ddG)</li>
<li>Binding affinity prediction: <code>FoldX</code>, <code>Rosetta interface_ddg</code>, AF2 PAE</li>
<li>Multiple "design-validate" iterations: experimental feedback returns to computation</li>
</ul>
<h2>6. Quick Reference Table for Tool Selection</h2>
<table>
<thead>
<tr>
<th>Task</th>
<th>Primary Tool</th>
<th>Alternative</th>
</tr>
</thead>
<tbody>
<tr>
<td>Structure prediction</td>
<td>AlphaFold2/3 (ColabFold)</td>
<td>ESMFold (speed)</td>
</tr>
<tr>
<td>Backbone generation</td>
<td>RFdiffusion</td>
<td>Chroma</td>
</tr>
<tr>
<td>Sequence design</td>
<td>ProteinMPNN</td>
<td>Rosetta FixBB</td>
</tr>
<tr>
<td>Thermostability design</td>
<td>Rosetta (ddg)</td>
<td>FoldX</td>
</tr>
<tr>
<td>Binding affinity assessment</td>
<td>Rosetta interface_ddg</td>
<td>FoldX / AF2 PAE</td>
</tr>
<tr>
<td>Interface design</td>
<td>Rosetta protocols</td>
<td>RFdiffusion + MPNN</td>
</tr>
<tr>
<td>Large-scale screening</td>
<td>ESMFold + ProteinMPNN</td>
<td>—</td>
</tr>
</tbody>
</table>
<h2>7. Experimental Validation Loop</h2>
<p>Computational design must return to experiments:</p>
<pre><code>Expression (E. coli / yeast) → Purification → Characterization
  → SEC / DSC (stability)
  → Surface plasmon resonance SPR / BLI (affinity)
  → Structural validation (crystallography / cryo-EM / AlphaFold comparison)
</code></pre>
<p><strong>Design success rate reference</strong>: Literature reports that AF2-guided RFdiffusion binder design can achieve experimental positive rates of 10–20% (far exceeding traditional methods), but every successful case involves multiple rounds of iteration.</p>
<h2>8. Summary</h2>
<ul>
<li><strong>Classic paradigm</strong>: Rosetta energy function (interpretable, adjustable)</li>
<li><strong>Prediction paradigm</strong>: AlphaFold family (validation, complex prediction)</li>
<li><strong>Generative paradigm</strong>: RFdiffusion + ProteinMPNN + AF2 validation loop (currently mainstream)</li>
<li>The ultimate measure of design success is always <strong>experimental validation</strong></li>
</ul>
<p>The next article will introduce mainstream tools for virtual screening: a complete toolchain for docking small molecules to targets.</p>`,Cb=`<h1>蛋白质设计主流工具</h1>
<p>蛋白质设计（Protein Design）是计算生物学的重要方向：给定目标功能，设计出能折叠成特定结构的氨基酸序列。本文梳理从经典能量函数到生成式 AI 的主流工具谱系。</p>
<h2>1. 设计问题的两种范式</h2>
<h3>1.1 序列设计（Sequence Design）</h3>
<p>给定目标结构（骨架），寻找最稳定的氨基酸序列：</p>
<pre><code>结构（backbone）→ 序列设计 → 氨基酸序列 → 实验验证
</code></pre>
<h3>1.2 结构设计（De novo Design）</h3>
<p>从零设计新结构/新功能：</p>
<pre><code>功能需求 → 结构生成 → 序列设计 → 实验验证
</code></pre>
<h2>2. Rosetta：能量函数范式</h2>
<p>Rosetta 是结构预测与设计领域的"老牌帝国"，基于物理与知识混合的能量函数。</p>
<h3>2.1 核心能量函数</h3>
<ul>
<li><strong>REF2015</strong>：经典全原子能量函数（范德华、氢键、溶剂化、静电、二面角偏好）</li>
<li>打分 = 各项能量加权和，越低越好</li>
<li>Rosetta 的"设计"本质是<strong>能量最小化</strong>：在序列/结构空间搜索低能态</li>
</ul>
<h3>2.2 常用协议</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 序列设计（FixBB：固定骨架设计）</span>
rosetta_scripts.default.linuxgccrelease \\
  -s input.pdb \\
  -parser:protocol fixbb.xml \\
  -out:prefix designed_

<span class="hljs-comment"># fixbb.xml 核心片段</span>
<span class="hljs-comment"># &#x3C;TaskOperations></span>
<span class="hljs-comment">#   &#x3C;RestrictToRepacking name="repack"/></span>
<span class="hljs-comment">#   &#x3C;OperateOnResidueSubset name="no_design_core" selector="..." /></span>
<span class="hljs-comment"># &#x3C;/TaskOperations></span>
<span class="hljs-comment"># &#x3C;RosettaScripts></span>
<span class="hljs-comment">#   &#x3C;PackRotamersMover name="pack" task_operations="repack"/></span>
<span class="hljs-comment"># &#x3C;/RosettaScripts></span>
</code></pre>
<h3>2.3 Rosetta 的经典应用</h3>
<ul>
<li><strong>酶设计</strong>：催化活性位点支架设计（如 Kemp 消除酶）</li>
<li><strong>蛋白-蛋白界面设计</strong>：结合特异性改造</li>
<li><strong>binder 设计</strong>：高亲和力结合蛋白</li>
<li><strong>热稳定性改造</strong>：计算突变提高 Tm</li>
</ul>
<h3>2.4 优点与局限</h3>
<table>
<thead>
<tr>
<th>优点</th>
<th>局限</th>
</tr>
</thead>
<tbody>
<tr>
<td>物理可解释、可精细调控</td>
<td>计算量大、需要经验调参</td>
</tr>
<tr>
<td>支持任意骨架与修饰</td>
<td>对柔性区域预测有限</td>
</tr>
<tr>
<td>生态成熟（文档、教程多）</td>
<td>学习曲线陡峭</td>
</tr>
</tbody>
</table>
<h2>3. AlphaFold：预测即设计辅助</h2>
<h3>3.1 AlphaFold2（2021）</h3>
<ul>
<li>端到端深度学习：序列 → 结构（MSA + 注意力机制）</li>
<li>精度：CASP14 原子级精度（GDT > 90）</li>
<li>对设计的价值：<strong>快速验证设计序列的可折叠性</strong></li>
</ul>
<h3>3.2 AlphaFold3（2024）</h3>
<ul>
<li>统一框架预测蛋白质、核酸、小分子、离子复合物</li>
<li>引入扩散模型与配对表示</li>
<li>意义：设计靶标（binder 与配体的复合物结构）可直接预测</li>
</ul>
<pre><code class="hljs language-bash"><span class="hljs-comment"># ColabFold：AlphaFold2 的轻量实现（GPU 友好）</span>
colabfold_batch input.fasta out_dir \\
  --num-recycle 3 --model-type alphafold2_multimer_v3
</code></pre>
<h3>3.3 关键用法：验证设计</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 设计序列的可信度指标</span>
<span class="hljs-comment"># pLDDT > 90：高置信</span>
<span class="hljs-comment"># PAE &#x3C; 5：结构域间相对位置可信</span>
<span class="hljs-comment"># 设计验证流程：序列 → AF2 预测 → 与目标结构比对（RMSD）</span>
</code></pre>
<h2>4. 生成式设计：RFdiffusion 与 ProteinMPNN</h2>
<h3>4.1 RFdiffusion（2023，Baker Lab）</h3>
<p>基于<strong>去噪扩散模型</strong>的结构生成器：</p>
<ul>
<li>输入：功能基序（motif）、形状约束、结合靶点</li>
<li>输出：全新蛋白质骨架（backbone）</li>
<li>应用：binder 设计、对称寡聚体、酶活性位点支架</li>
</ul>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RFdiffusion 示例：设计结合蛋白</span>
run_inference.py \\
  scaffoldguided.scaffoldguided=True \\
  scaff_loader.contigs_map.json=contigs.json \\
  inference.output_prefix=outputs/design \\
  inference.num_designs=100 \\
  potentials.guidance_scale=3.0
</code></pre>
<h3>4.2 ProteinMPNN（2022）</h3>
<p><strong>序列设计神经网络</strong>：给定骨架 → 序列，与 RFdiffusion 形成"生成-编码"闭环：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># ProteinMPNN 序列设计</span>
python protein_mpnn_run.py \\
  --pdbpath scaffold.pdb \\
  --out_folder mpnn_out \\
  --num_seq_per_target 8 \\
  --batch_size 4
</code></pre>
<p><strong>标准工作流（Baker Lab 范式）</strong>：</p>
<pre><code>① RFdiffusion：生成骨架
   ↓
② ProteinMPNN：骨架 → 多组候选序列（采样温度 0.1-0.3）
   ↓
③ AF2 反向验证：序列 → 预测结构 → 与骨架 RMSD 比对
   ↓
④ 筛选高 pLDDT / 低 RMSD 候选
   ↓
⑤ 实验表达验证（酵母/噬菌体展示、ELISA）
</code></pre>
<h3>4.3 其他生成式工具</h3>
<ul>
<li><strong>Chroma</strong>：另一个扩散模型（生成 + 序列设计一体）</li>
<li><strong>Genie</strong>：更快的蛋白质扩散模型</li>
<li><strong>ESMFold</strong>（Meta，2022）：无 MSA 的极速预测（每序列 &#x3C; 1 秒），适合大规模设计筛选</li>
<li><strong>ESM3</strong>：多模态生成（序列+结构+功能）</li>
</ul>
<h2>5. binder 设计专题</h2>
<p>Binder（结合蛋白）设计是当前最活跃的应用场景：</p>
<h3>5.1 经典流程（结合现有工具）</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. 靶点表面热点分析（Rosetta）</span>
<span class="hljs-comment"># 2. 热区引导的 binder 骨架生成（RFdiffusion / 或 hotspot 约束）</span>
<span class="hljs-comment"># 3. ProteinMPNN 序列设计</span>
<span class="hljs-comment"># 4. AF2 复合物预测验证（binder + target 对接）</span>
<span class="hljs-comment"># 5. 实验筛选（酵母展示、FACS 分选）</span>
</code></pre>
<h3>5.2 关键考量</h3>
<ul>
<li>结合热点残基（hotspot）的确定：实验突变扫描（alanine scan）或计算（Rosetta ddG）</li>
<li>结合亲和力预测：<code>FoldX</code>、<code>Rosetta interface_ddg</code>、AF2 PAE</li>
<li>多轮"设计-验证"迭代：实验反馈回到计算</li>
</ul>
<h2>6. 工具选择速查表</h2>
<table>
<thead>
<tr>
<th>任务</th>
<th>首选工具</th>
<th>备选</th>
</tr>
</thead>
<tbody>
<tr>
<td>结构预测</td>
<td>AlphaFold2/3（ColabFold）</td>
<td>ESMFold（速度）</td>
</tr>
<tr>
<td>骨架生成</td>
<td>RFdiffusion</td>
<td>Chroma</td>
</tr>
<tr>
<td>序列设计</td>
<td>ProteinMPNN</td>
<td>Rosetta FixBB</td>
</tr>
<tr>
<td>热稳定性设计</td>
<td>Rosetta（ddg）</td>
<td>FoldX</td>
</tr>
<tr>
<td>结合亲和力评估</td>
<td>Rosetta interface_ddg</td>
<td>FoldX / AF2 PAE</td>
</tr>
<tr>
<td>界面设计</td>
<td>Rosetta 协议</td>
<td>RFdiffusion + MPNN</td>
</tr>
<tr>
<td>大规模筛选</td>
<td>ESMFold + ProteinMPNN</td>
<td>—</td>
</tr>
</tbody>
</table>
<h2>7. 实验验证闭环</h2>
<p>计算设计必须回到实验：</p>
<pre><code>表达（E. coli / 酵母）→ 纯化 → 表征
  → SEC / DSC（稳定性）
  → 表面等离子共振 SPR / BLI（亲和力）
  → 结构验证（晶体 / 冷冻电镜 / AlphaFold 对照）
</code></pre>
<p><strong>设计成功率参考</strong>：文献报道 AF2 引导的 RFdiffusion binder 设计，实验阳性率可达 10–20%（远超传统方法），但每个成功案例背后都有多轮迭代。</p>
<h2>8. 小结</h2>
<ul>
<li><strong>经典范式</strong>：Rosetta 能量函数（可解释、可调控）</li>
<li><strong>预测范式</strong>：AlphaFold 系列（验证、复合物预测）</li>
<li><strong>生成范式</strong>：RFdiffusion + ProteinMPNN + AF2 验证闭环（当前主流）</li>
<li>设计成功的衡量标准始终是<strong>实验验证</strong></li>
</ul>
<p>下一篇将介绍虚拟筛选的主流工具：把小分子与靶点对接的完整工具链。</p>`,kb=`<h1>Mainstream Tools for Virtual Screening</h1>
<p>Virtual screening (VS) uses computation to search for potential active molecules in million-scale compound libraries and is a core approach in the early stage of drug discovery. This article introduces mainstream docking software, the complete workflow, and important considerations.</p>
<h2>1. Basic Logic of Virtual Screening</h2>
<pre><code>Receptor structure (protein/nucleic acid)
   ↓
① Receptor preparation (add hydrogens, protonation, grid)
   ↓
② Compound library (ZINC / Enamine / self-built library)
   ↓
③ Ligand preparation (3D conformation, protonation, force field parameters)
   ↓
④ Molecular docking (conformational sampling + scoring)
   ↓
⑤ Ranking and selection (Top 100-1000)
   ↓
⑥ Experimental validation (IC50 / activity assay)
</code></pre>
<h2>2. Mainstream Docking Software</h2>
<h3>2.1 AutoDock Vina (Top Open-Source Choice)</h3>
<p>The most popular open-source docking software, offering a good balance between speed and accuracy:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. Receptor and ligand preparation (using MGLTools or scripts)</span>
<span class="hljs-comment">#    Receptor: pdb → pdbqt (add hydrogens, merge nonpolar hydrogens)</span>
<span class="hljs-comment">#    Ligand: sdf/mol2 → pdbqt</span>

<span class="hljs-comment"># 2. Configuration file (conf.txt)</span>
receptor = receptor.pdbqt
ligand = ligand.pdbqt
center_x = 10.5
center_y = 20.3
center_z = -5.8
size_x = 22
size_y = 22
size_z = 22
exhaustiveness = 16
num_modes = 9

<span class="hljs-comment"># 3. Run</span>
vina --config conf.txt --out results.pdbqt --<span class="hljs-built_in">log</span> log.txt

<span class="hljs-comment"># 4. Interpretation: the lower the affinity (kcal/mol), the better</span>
<span class="hljs-comment">#  mode |   affinity | dist from best mode</span>
<span class="hljs-comment"># -----+------------+----------------------</span>
<span class="hljs-comment">#    1 |   -9.2     | 0.000</span>
<span class="hljs-comment">#    2 |   -8.7     | 2.315</span>
</code></pre>
<h3>2.2 Glide (Schrödinger, Commercial)</h3>
<p>The standard in the pharmaceutical industry:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Prepare workflow (Maestro GUI or command line)</span>
glide -overwrite -adjust \\
  -receptor receptor.maegz \\
  -ligand ligands.maegz \\
  -p glide_docking.inp \\
  -NJOBS 16

<span class="hljs-comment"># Three precision modes</span>
<span class="hljs-comment"># HTVS (high-throughput, fast) → SP (standard precision) → XP (high precision, slow)</span>
</code></pre>
<p><strong>Advantages</strong>: Accurate flexible ligand handling, well-calibrated scoring, and support from the large pharma ecosystem.</p>
<h3>2.3 DOCK 6 (Open Source)</h3>
<p>Classic geometric matching algorithm:</p>
<pre><code class="hljs language-bash">dock6 -i dock.in -o dock.out
<span class="hljs-comment"># Sphere generation → orientation search → scoring</span>
<span class="hljs-comment"># Strengths: binding pocket shape matching, initial screening of large libraries</span>
</code></pre>
<h3>2.4 Other Commonly Used Software</h3>
<table>
<thead>
<tr>
<th>Software</th>
<th>Features</th>
</tr>
</thead>
<tbody>
<tr>
<td><strong>AutoDock4</strong></td>
<td>Classic genetic algorithm, high accuracy but slow</td>
</tr>
<tr>
<td><strong>rDock</strong></td>
<td>Open source, suitable for large libraries, easy to script</td>
</tr>
<tr>
<td><strong>GOLD</strong></td>
<td>Genetic algorithm, developed by CCDC</td>
</tr>
<tr>
<td><strong>LeDock</strong></td>
<td>Fast and easy to use (free for academia)</td>
</tr>
<tr>
<td><strong>Plants</strong></td>
<td>Based on Ant Colony Optimization</td>
</tr>
<tr>
<td><strong>smina</strong></td>
<td>Enhanced version of Vina (supports custom scoring)</td>
</tr>
</tbody>
</table>
<h2>3. Receptor Preparation (Critical First Step)</h2>
<h3>3.1 Structure Sources</h3>
<ul>
<li>Crystal/cryo-EM structures (PDB)</li>
<li>AlphaFold-predicted structures (when no experimental structure is available)</li>
<li>Note: handle flexibility (side chains) near the binding pocket</li>
</ul>
<h3>3.2 Key Preparation Points</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Common tools</span>
<span class="hljs-comment"># Schrödinger Protein Prep Wizard (commercial)</span>
<span class="hljs-comment"># ADFRsuite / MGLTools prepare_receptor4.py (free)</span>

python prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt

<span class="hljs-comment"># Key steps</span>
<span class="hljs-comment"># 1. Remove water molecules (except conserved waters)</span>
<span class="hljs-comment"># 2. Add hydrogens, determine protonation states (around pH 7.4)</span>
<span class="hljs-comment"># 3. Assign bond orders and charges (Amber/CHARMM force fields)</span>
<span class="hljs-comment"># 4. Define the binding pocket (known active site or cavity detection: fpocket / DoGSiteScorer)</span>
</code></pre>
<h2>4. Ligand Libraries and Preparation</h2>
<h3>4.1 Compound Library Sources</h3>
<table>
<thead>
<tr>
<th>Library</th>
<th>Scale</th>
<th>Features</th>
</tr>
</thead>
<tbody>
<tr>
<td>ZINC20</td>
<td>Billions</td>
<td>Free, downloadable, with 3D conformations</td>
</tr>
<tr>
<td>Enamine REAL</td>
<td>> 40 billion</td>
<td>Highly synthesizable, commercially available</td>
</tr>
<tr>
<td>ChEMBL</td>
<td>Bioactivity data</td>
<td>Known active compounds</td>
</tr>
<tr>
<td>Self-built library</td>
<td>Custom</td>
<td>Derivative design</td>
</tr>
</tbody>
</table>
<h3>4.2 Ligand Preparation</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 3D conformation generation</span>
obabel ligand.smi -O ligand.sdf --gen3d -p 7.4

<span class="hljs-comment"># Multiple conformations (flexibility is important, especially for ring systems)</span>
obabel ligand.sdf -O conf.sdf --conformer --nconf 50

<span class="hljs-comment"># Protonation and charges</span>
<span class="hljs-comment"># Schrödinger LigPrep / Open Babel / RDKit</span>

<span class="hljs-comment"># RDKit preparation (Python)</span>
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

mol = Chem.MolFromSmiles(<span class="hljs-string">"CC(=O)Oc1ccccc1C(=O)O"</span>)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
</code></pre>
<h2>5. Scoring Functions</h2>
<h3>5.1 Three Types of Scoring Functions</h3>
<table>
<thead>
<tr>
<th>Type</th>
<th>Representatives</th>
<th>Principle</th>
</tr>
</thead>
<tbody>
<tr>
<td>Force-field-based</td>
<td>DOCK, AutoDock4</td>
<td>van der Waals + electrostatics + internal energy</td>
</tr>
<tr>
<td>Empirical</td>
<td>Glide, Vina</td>
<td>Statistical fitting to experimental data</td>
</tr>
<tr>
<td>Knowledge-based</td>
<td>GOLD/ASP</td>
<td>Statistics of atom-pair distance distributions</td>
</tr>
</tbody>
</table>
<h3>5.2 Limitations of Scoring Functions (Must-Know)</h3>
<ul>
<li><strong>Not precise</strong>: correlation coefficients are typically 0.3–0.6; they can only rank, not provide absolute affinities</li>
<li>Entropy effects, solvent effects, and induced fit are difficult to describe accurately</li>
<li><strong>Recommendation</strong>: docking scores + cross-validation with multiple tools + refined calculations using molecular dynamics (MM/PBSA, FEP)</li>
</ul>
<h2>6. Complete High-Throughput Screening Workflow (Hands-On HTS)</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. Library preparation: split SDF + convert to pdbqt</span>
<span class="hljs-built_in">mkdir</span> -p ligands_pdbqt
obabel library.sdf -O ligands_pdbqt/lig_.pdbqt -m -p 7.4

<span class="hljs-comment"># 2. Batch docking (smina/Vina multi-file mode)</span>
smina -r receptor.pdbqt -l library.sdf \\
  --center_x 10.5 --center_y 20.3 --center_z -5.8 \\
  --size_x 22 --size_y 22 --size_z 22 \\
  --num_modes 1 --cpu 32 \\
  -o docked.sdf --<span class="hljs-built_in">log</span> scores.txt

<span class="hljs-comment"># 3. Ranking (by affinity)</span>
<span class="hljs-built_in">sort</span> -k2 -n scores.txt | <span class="hljs-built_in">head</span> -100
</code></pre>
<h3>6.1 Multi-Stage Funnel Strategy</h3>
<pre><code>Coarse screening (HTVS / Vina, million-scale) → top 5-10%
  ↓
Standard precision (SP / multiple conformations) → top 10-20%
  ↓
High precision (XP / induced fit) → Top 100-1000
  ↓
Consensus scoring (cross-validation with multiple software tools) → Top 50-200
  ↓
Visual inspection (reasonableness of binding modes) → Top 10-50
  ↓
Experimental testing
</code></pre>
<h3>6.2 Consensus Scoring</h3>
<p>Scoring functions from different software tools complement each other; taking intersections/ranking averages can significantly improve enrichment rates:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Example: take the intersection of Top 100 from Vina + Glide + LeDock</span>
<span class="hljs-comment"># Or average after normalizing ranks</span>
</code></pre>
<h2>7. Validation and Advanced Methods</h2>
<h3>7.1 Benchmarking</h3>
<ul>
<li><strong>DUD-E</strong>: a standard dataset (containing actives and decoys) used to evaluate a software's ability to discriminate</li>
<li>Metrics: AUC, enrichment factors (EF1%, EF5%)</li>
</ul>
<h3>7.2 Molecular Dynamics Post-Processing</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Docking results → MD simulation (GROMACS/AMBER) → MM/PBSA binding free energy</span>
gmx mdrun -deffnm complex
gmx_MMPBSA -O -i mmpbsa.in -cs complex.tpr -ct traj.xtc ...
</code></pre>
<h3>7.3 Combining Structure-Based and Ligand-Based Approaches</h3>
<ul>
<li><strong>Structure-based</strong> (SBDD): docking + MD + FEP</li>
<li><strong>Ligand-based</strong> (LBDD): pharmacophore, QSAR, similarity search</li>
<li>Combination strategies are the mainstream of modern virtual screening</li>
</ul>
<h2>8. Common Pitfalls</h2>
<table>
<thead>
<tr>
<th>Pitfall</th>
<th>Consequence</th>
<th>Avoidance</th>
</tr>
</thead>
<tbody>
<tr>
<td>Incorrect pocket definition</td>
<td>Entire library misaligned</td>
<td>Use a known active site or experimental information</td>
</tr>
<tr>
<td>Rigid receptor</td>
<td>Misses induced fit</td>
<td>Flexible side chains / ensemble docking</td>
</tr>
<tr>
<td>Incorrect protonation states</td>
<td>Distorted electrostatic scoring</td>
<td>pH-dependent preparation (pH 7.4)</td>
</tr>
<tr>
<td>Overconfidence in scoring functions</td>
<td>High false-positive rate</td>
<td>Consensus scoring + visual inspection</td>
</tr>
<tr>
<td>Ignoring synthesizability</td>
<td>Hits cannot be synthesized</td>
<td>Filter PAINS and drug-likeness rules (Lipinski/Veber)</td>
</tr>
</tbody>
</table>
<h2>9. Summary</h2>
<ul>
<li>Three tiers: <strong>Vina</strong> (open source, fast) → <strong>Glide</strong> (commercial standard) → <strong>DOCK/rDock</strong> (large libraries)</li>
<li>Workflow: receptor preparation → ligand library → docking → multi-stage funnel → experimental validation</li>
<li>Scoring functions can only rank; consensus scoring + MD refinement improve reliability</li>
<li>The end point of virtual screening is always experimentation</li>
</ul>
<p>This completes the computational biology direction: protein design tools + virtual screening tools, with two main threads forming a closed loop of "design-screen" computational drug discovery.</p>`,xb=`<h1>虚拟筛选主流工具</h1>
<p>虚拟筛选（Virtual Screening, VS）通过计算在百万级化合物库中寻找潜在活性分子，是药物发现早期阶段的核心手段。本文介绍主流对接软件、完整流程与注意事项。</p>
<h2>1. 虚拟筛选的基本逻辑</h2>
<pre><code>受体结构（蛋白/核酸）
   ↓
① 受体准备（加氢、质子化、网格）
   ↓
② 化合物库（ZINC / Enamine / 自建库）
   ↓
③ 配体准备（3D 构象、质子化、力场参数）
   ↓
④ 分子对接（构象采样 + 打分）
   ↓
⑤ 排序筛选（Top 100-1000）
   ↓
⑥ 实验验证（IC50 / 活性测定）
</code></pre>
<h2>2. 主流对接软件</h2>
<h3>2.1 AutoDock Vina（开源首选）</h3>
<p>最流行的开源对接软件，速度与精度平衡良好：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. 受体与配体准备（用 MGLTools 或脚本）</span>
<span class="hljs-comment">#    受体：pdb → pdbqt（加氢、合并非极性氢）</span>
<span class="hljs-comment">#    配体：sdf/mol2 → pdbqt</span>

<span class="hljs-comment"># 2. 配置文件（conf.txt）</span>
receptor = receptor.pdbqt
ligand = ligand.pdbqt
center_x = 10.5
center_y = 20.3
center_z = -5.8
size_x = 22
size_y = 22
size_z = 22
exhaustiveness = 16
num_modes = 9

<span class="hljs-comment"># 3. 运行</span>
vina --config conf.txt --out results.pdbqt --<span class="hljs-built_in">log</span> log.txt

<span class="hljs-comment"># 4. 解读：affinity（kcal/mol）越低越好</span>
<span class="hljs-comment">#  mode |   affinity | dist from best mode</span>
<span class="hljs-comment"># -----+------------+----------------------</span>
<span class="hljs-comment">#    1 |   -9.2     | 0.000</span>
<span class="hljs-comment">#    2 |   -8.7     | 2.315</span>
</code></pre>
<h3>2.2 Glide（Schrödinger，商业）</h3>
<p>制药行业标准：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 准备工作流（Maestro GUI 或命令行）</span>
glide -overwrite -adjust \\
  -receptor receptor.maegz \\
  -ligand ligands.maegz \\
  -p glide_docking.inp \\
  -NJOBS 16

<span class="hljs-comment"># 三种精度模式</span>
<span class="hljs-comment"># HTVS（高通量，快）→ SP（标准精度）→ XP（高精度，慢）</span>
</code></pre>
<p><strong>优势</strong>：精确的配体柔性处理、完善的打分校准、大型药企生态支持。</p>
<h3>2.3 DOCK 6（开源）</h3>
<p>经典几何匹配算法（shape matching）：</p>
<pre><code class="hljs language-bash">dock6 -i dock.in -o dock.out
<span class="hljs-comment"># 球集（sphere）生成 → 取向搜索 → 打分</span>
<span class="hljs-comment"># 擅长：结合口袋形状匹配、大库初筛</span>
</code></pre>
<h3>2.4 其他常用软件</h3>
<table>
<thead>
<tr>
<th>软件</th>
<th>特点</th>
</tr>
</thead>
<tbody>
<tr>
<td><strong>AutoDock4</strong></td>
<td>经典遗传算法，精度高速度慢</td>
</tr>
<tr>
<td><strong>rDock</strong></td>
<td>开源、适合大库、易脚本化</td>
</tr>
<tr>
<td><strong>GOLD</strong></td>
<td>遗传算法，CCDC 出品</td>
</tr>
<tr>
<td><strong>LeDock</strong></td>
<td>快速、易用（学术免费）</td>
</tr>
<tr>
<td><strong>Plants</strong></td>
<td>基于 Ant Colony 优化</td>
</tr>
<tr>
<td><strong>smina</strong></td>
<td>Vina 的增强版（支持自定义打分）</td>
</tr>
</tbody>
</table>
<h2>3. 受体准备（关键第一步）</h2>
<h3>3.1 结构来源</h3>
<ul>
<li>晶体/冷冻电镜结构（PDB）</li>
<li>AlphaFold 预测结构（无实验结构时）</li>
<li>注意：结合口袋附近的柔性（侧链）处理</li>
</ul>
<h3>3.2 准备要点</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 常用工具</span>
<span class="hljs-comment"># Schrödinger Protein Prep Wizard（商业）</span>
<span class="hljs-comment"># ADFRsuite / MGLTools prepare_receptor4.py（免费）</span>

python prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt

<span class="hljs-comment"># 关键步骤</span>
<span class="hljs-comment"># 1. 去除水分子（保守水除外）</span>
<span class="hljs-comment"># 2. 加氢、确定质子化状态（pH 7.4 附近）</span>
<span class="hljs-comment"># 3. 分配键序与电荷（Amber/CHARMM 力场）</span>
<span class="hljs-comment"># 4. 定义结合口袋（已知活性位点 or 空腔检测：fpocket / DoGSiteScorer）</span>
</code></pre>
<h2>4. 配体库与准备</h2>
<h3>4.1 化合物库来源</h3>
<table>
<thead>
<tr>
<th>库</th>
<th>规模</th>
<th>特点</th>
</tr>
</thead>
<tbody>
<tr>
<td>ZINC20</td>
<td>数十亿</td>
<td>免费、可下载、带 3D 构象</td>
</tr>
<tr>
<td>Enamine REAL</td>
<td>> 400 亿</td>
<td>可合成性高、商购</td>
</tr>
<tr>
<td>ChEMBL</td>
<td>活性数据</td>
<td>已知活性化合物</td>
</tr>
<tr>
<td>自建库</td>
<td>定制</td>
<td>衍生物设计</td>
</tr>
</tbody>
</table>
<h3>4.2 配体准备</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 3D 构象生成</span>
obabel ligand.smi -O ligand.sdf --gen3d -p 7.4

<span class="hljs-comment"># 多构象（柔性重要，尤其环状体系）</span>
obabel ligand.sdf -O conf.sdf --conformer --nconf 50

<span class="hljs-comment"># 质子化与电荷</span>
<span class="hljs-comment"># Schrödinger LigPrep / Open Babel / RDKit</span>

<span class="hljs-comment"># RDKit 准备（Python）</span>
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

mol = Chem.MolFromSmiles(<span class="hljs-string">"CC(=O)Oc1ccccc1C(=O)O"</span>)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
</code></pre>
<h2>5. 打分函数（Scoring Functions）</h2>
<h3>5.1 三类打分函数</h3>
<table>
<thead>
<tr>
<th>类型</th>
<th>代表</th>
<th>原理</th>
</tr>
</thead>
<tbody>
<tr>
<td>力场型</td>
<td>DOCK、AutoDock4</td>
<td>范德华+静电+内能</td>
</tr>
<tr>
<td>经验型</td>
<td>Glide、Vina</td>
<td>统计拟合实验数据</td>
</tr>
<tr>
<td>知识型</td>
<td>GOLD/ASP</td>
<td>原子对距离分布统计</td>
</tr>
</tbody>
</table>
<h3>5.2 打分函数的局限（必须知道）</h3>
<ul>
<li><strong>不精确</strong>：相关系数通常 0.3–0.6，只能排序不能给绝对亲和力</li>
<li>熵效应、溶剂效应、诱导拟合难以准确描述</li>
<li><strong>建议</strong>：对接分数 + 多工具交叉验证 + 分子动力学（MM/PBSA、FEP）精算</li>
</ul>
<h2>6. 完整高通量筛选流程（HTS 实操）</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. 库准备：SDF 拆分 + 转 pdbqt</span>
<span class="hljs-built_in">mkdir</span> -p ligands_pdbqt
obabel library.sdf -O ligands_pdbqt/lig_.pdbqt -m -p 7.4

<span class="hljs-comment"># 2. 批量对接（smina/Vina 多文件模式）</span>
smina -r receptor.pdbqt -l library.sdf \\
  --center_x 10.5 --center_y 20.3 --center_z -5.8 \\
  --size_x 22 --size_y 22 --size_z 22 \\
  --num_modes 1 --cpu 32 \\
  -o docked.sdf --<span class="hljs-built_in">log</span> scores.txt

<span class="hljs-comment"># 3. 排序（按 affinity）</span>
<span class="hljs-built_in">sort</span> -k2 -n scores.txt | <span class="hljs-built_in">head</span> -100
</code></pre>
<h3>6.1 多阶段漏斗策略</h3>
<pre><code>粗筛（HTVS / Vina，百万级）→ 前 5-10%
  ↓
标准精度（SP / 多构象）→ 前 10-20%
  ↓
高精度（XP / 诱导拟合）→ Top 100-1000
  ↓
共识打分（多个软件交叉）→ Top 50-200
  ↓
视觉检查（结合模式合理性）→ Top 10-50
  ↓
实验测试
</code></pre>
<h3>6.2 共识打分（Consensus Scoring）</h3>
<p>不同软件打分函数互补，取交集/排名平均可显著提升富集率：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 示例：Vina + Glide + LeDock 三软件 Top 100 取交集</span>
<span class="hljs-comment"># 或按排名归一化后取平均</span>
</code></pre>
<h2>7. 验证与高级方法</h2>
<h3>7.1 基准测试（Benchmark）</h3>
<ul>
<li><strong>DUD-E</strong>：标准数据集（含活性物与诱饵 decoys），评估软件区分能力</li>
<li>指标：AUC、富集因子（EF1%、EF5%）</li>
</ul>
<h3>7.2 分子动力学后处理</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 对接结果 → MD 模拟（GROMACS/AMBER）→ MM/PBSA 结合自由能</span>
gmx mdrun -deffnm complex
gmx_MMPBSA -O -i mmpbsa.in -cs complex.tpr -ct traj.xtc ...
</code></pre>
<h3>7.3 结构基与配体基结合</h3>
<ul>
<li><strong>基于结构</strong>（SBDD）：对接 + MD + FEP</li>
<li><strong>基于配体</strong>（LBDD）：药效团、QSAR、相似性搜索</li>
<li>组合策略是现代虚拟筛选主流</li>
</ul>
<h2>8. 常见陷阱</h2>
<table>
<thead>
<tr>
<th>陷阱</th>
<th>后果</th>
<th>规避</th>
</tr>
</thead>
<tbody>
<tr>
<td>口袋定义错误</td>
<td>全库错位</td>
<td>用已知活性位点或实验信息</td>
</tr>
<tr>
<td>受体刚性</td>
<td>错过诱导拟合</td>
<td>柔性侧链 / 集合对接（ensemble docking）</td>
</tr>
<tr>
<td>质子化状态错误</td>
<td>静电打分失真</td>
<td>pH 相关准备（pH 7.4）</td>
</tr>
<tr>
<td>打分函数过信</td>
<td>假阳性率高</td>
<td>共识打分 + 视觉检查</td>
</tr>
<tr>
<td>忽略可合成性</td>
<td>筛出无法合成</td>
<td>过滤 PAINS、类药性规则（Lipinski/Veber）</td>
</tr>
</tbody>
</table>
<h2>9. 小结</h2>
<ul>
<li>三梯队：<strong>Vina</strong>（开源快速）→ <strong>Glide</strong>（商业标准）→ <strong>DOCK/rDock</strong>（大库）</li>
<li>流程：受体准备 → 配体库 → 对接 → 多阶段漏斗 → 实验验证</li>
<li>打分函数只能排序；共识打分 + MD 精算提升可信度</li>
<li>虚拟筛选的终点永远是实验</li>
</ul>
<p>至此计算生物学方向完成：蛋白质设计工具 + 虚拟筛选工具，两条主线构成"设计-筛选"的计算药物发现闭环。</p>`,Pb=`<h1>Bash Programming: Variables, Conditionals, Loops, and Practical Scripts</h1>
<p>Bash scripts organize command-line operations into reusable programs and are a core tool for automated batch processing. This article covers variables, conditionals, loops, functions, and common tasks, from script basics to practical applications.</p>
<h2>1. Script Basics</h2>
<h3>1.1 Script Structure</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># The first line shebang: specifies the interpreter</span>

<span class="hljs-comment"># Give the script execute permission</span>
<span class="hljs-built_in">chmod</span> +x myscript.sh
<span class="hljs-comment"># Run</span>
./myscript.sh
<span class="hljs-comment"># Or run directly with bash (no execute permission needed)</span>
bash myscript.sh
</code></pre>
<h3>1.2 Debugging</h3>
<pre><code class="hljs language-bash">bash -x script.sh    <span class="hljs-comment"># Trace execution, display each command</span>
<span class="hljs-comment"># Enable inside the script: set -x (trace) / set -e (exit on error) / set -u (error on undefined variable)</span>
</code></pre>
<p>For production scripts, it's recommended to add at the beginning:</p>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-built_in">set</span> -euo pipefail   <span class="hljs-comment"># Exit on error, error on undefined variables, pipe failure detection</span>
</code></pre>
<h2>2. Variables</h2>
<h3>2.1 Definition and Reference</h3>
<pre><code class="hljs language-bash">name=<span class="hljs-string">"zorrooz"</span>        <span class="hljs-comment"># Note: no spaces around =</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Hello, <span class="hljs-variable">$name</span>"</span>   <span class="hljs-comment"># Variable expansion inside double quotes</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">'Hello, $name'</span>   <span class="hljs-comment"># Literal output inside single quotes</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Length: <span class="hljs-variable">\${#name}</span>"</span> <span class="hljs-comment"># String length</span>
</code></pre>
<h3>2.2 Command Substitution and Arithmetic</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Command substitution: $(...)</span>
today=$(<span class="hljs-built_in">date</span> +%F)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Today is <span class="hljs-variable">$today</span>"</span>
files=$(<span class="hljs-built_in">ls</span> | <span class="hljs-built_in">wc</span> -l)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Total: <span class="hljs-variable">$files</span> files"</span>

<span class="hljs-comment"># Arithmetic: $((...))</span>
a=10
b=3
<span class="hljs-built_in">echo</span> $((a + b))    <span class="hljs-comment"># 13</span>
<span class="hljs-built_in">echo</span> $((a % b))    <span class="hljs-comment"># 1</span>
</code></pre>
<h3>2.3 Positional Parameters</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Script name: <span class="hljs-variable">$0</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"First argument: <span class="hljs-variable">$1</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Second argument: <span class="hljs-variable">$2</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"All arguments: <span class="hljs-variable">$@</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Argument count: <span class="hljs-variable">$#</span>"</span>

<span class="hljs-comment"># Run: ./script.sh file1 file2</span>
</code></pre>
<h3>2.4 Special Variables</h3>
<pre><code class="hljs language-bash">$?    <span class="hljs-comment"># Exit code of the last command (0 success, non-zero failure)</span>
$$    <span class="hljs-comment"># PID of the current process</span>
$!    <span class="hljs-comment"># PID of the background process</span>
</code></pre>
<h2>3. Conditional Statements</h2>
<h3>3.1 if / elif / else</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
score=<span class="hljs-variable">$1</span>

<span class="hljs-keyword">if</span> (( score >= <span class="hljs-number">90</span> )); <span class="hljs-keyword">then</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Excellent"</span>
<span class="hljs-keyword">elif</span> (( score >= <span class="hljs-number">60</span> )); <span class="hljs-keyword">then</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Pass"</span>
<span class="hljs-keyword">else</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Fail"</span>
<span class="hljs-keyword">fi</span>
</code></pre>
<h3>3.2 test Expressions</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># File tests</span>
[ -f file.txt ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Is a regular file"</span>
[ -d <span class="hljs-built_in">dir</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Is a directory"</span>
[ -e path ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Exists"</span>
[ -s file ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Is a non-empty file"</span>

<span class="hljs-comment"># String comparison</span>
[ <span class="hljs-string">"<span class="hljs-variable">$name</span>"</span> = <span class="hljs-string">"zorrooz"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Match"</span>
[ -z <span class="hljs-string">"<span class="hljs-variable">$var</span>"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Variable is empty"</span>     <span class="hljs-comment"># -n non-empty</span>

<span class="hljs-comment"># Numeric comparison (inside [ ])</span>
[ <span class="hljs-string">"<span class="hljs-variable">$a</span>"</span> -gt <span class="hljs-string">"<span class="hljs-variable">$b</span>"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"a > b"</span>    <span class="hljs-comment"># -eq -ne -gt -ge -lt -le</span>

<span class="hljs-comment"># Logical combination</span>
[ -f <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] &#x26;&#x26; [ -s <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Non-empty file"</span>
[ -f <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] || <span class="hljs-built_in">echo</span> <span class="hljs-string">"File does not exist"</span>
</code></pre>
<h3>3.3 case Branches</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-keyword">case</span> <span class="hljs-string">"<span class="hljs-variable">$1</span>"</span> <span class="hljs-keyword">in</span>
  start)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Start the service"</span>
    ;;
  stop|<span class="hljs-built_in">kill</span>)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Stop the service"</span>
    ;;
  *)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Usage: <span class="hljs-variable">$0</span> {start|stop}"</span>
    <span class="hljs-built_in">exit</span> 1
    ;;
<span class="hljs-keyword">esac</span>
</code></pre>
<h2>4. Loops</h2>
<h3>4.1 for Loops</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Iterate over a list</span>
<span class="hljs-keyword">for</span> fruit <span class="hljs-keyword">in</span> apple banana cherry; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Processing <span class="hljs-variable">$fruit</span>"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-comment"># Iterate over a range</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> {1..5}; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Iteration <span class="hljs-variable">$i</span>"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-comment"># Iterate over files (most common)</span>
<span class="hljs-keyword">for</span> file <span class="hljs-keyword">in</span> *.fasta; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Processing <span class="hljs-variable">$file</span>"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-comment"># C-style</span>
<span class="hljs-keyword">for</span> ((i = <span class="hljs-number">0</span>; i &#x3C; <span class="hljs-number">10</span>; i++)); <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-variable">$i</span>
<span class="hljs-keyword">done</span>
</code></pre>
<h3>4.2 while Loops</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Read a file line by line</span>
<span class="hljs-keyword">while</span> IFS= <span class="hljs-built_in">read</span> -r line; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Read: <span class="hljs-variable">$line</span>"</span>
<span class="hljs-keyword">done</span> &#x3C; input.txt

<span class="hljs-comment"># Counting loop</span>
count=0
<span class="hljs-keyword">while</span> [ <span class="hljs-variable">$count</span> -lt 5 ]; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-variable">$count</span>
    count=$((count + <span class="hljs-number">1</span>))
<span class="hljs-keyword">done</span>
</code></pre>
<h3>4.3 break / continue</h3>
<pre><code class="hljs language-bash"><span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> {1..10}; <span class="hljs-keyword">do</span>
    [ <span class="hljs-variable">$i</span> -eq 3 ] &#x26;&#x26; <span class="hljs-built_in">continue</span>   <span class="hljs-comment"># Skip 3</span>
    [ <span class="hljs-variable">$i</span> -eq 7 ] &#x26;&#x26; <span class="hljs-built_in">break</span>      <span class="hljs-comment"># Terminate at 7</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-variable">$i</span>
<span class="hljs-keyword">done</span>
</code></pre>
<h2>5. Functions</h2>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>

<span class="hljs-comment"># Define a function</span>
<span class="hljs-function"><span class="hljs-title">greet</span></span>() {
    <span class="hljs-built_in">local</span> name=<span class="hljs-variable">$1</span>           <span class="hljs-comment"># local: local variable</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Hello, <span class="hljs-variable">$name</span>!"</span>
    <span class="hljs-built_in">return</span> 0                <span class="hljs-comment"># Return value (exit code)</span>
}

<span class="hljs-comment"># Call</span>
greet <span class="hljs-string">"zorrooz"</span>

<span class="hljs-comment"># Function returning a value (note: not return, but echo capture)</span>
<span class="hljs-function"><span class="hljs-title">add</span></span>() {
    <span class="hljs-built_in">echo</span> $(( <span class="hljs-variable">$1</span> + <span class="hljs-variable">$2</span> ))
}
result=$(add 3 4)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"3 + 4 = <span class="hljs-variable">$result</span>"</span>
</code></pre>
<h2>6. Arrays and Strings</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Arrays</span>
genes=(<span class="hljs-string">"TP53"</span> <span class="hljs-string">"BRCA1"</span> <span class="hljs-string">"EGFR"</span>)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${genes[0]}</span>"</span>          <span class="hljs-comment"># TP53</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${genes[@]}</span>"</span>          <span class="hljs-comment"># All elements</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${#genes[@]}</span>"</span>         <span class="hljs-comment"># Length</span>
genes+=(<span class="hljs-string">"MYC"</span>)              <span class="hljs-comment"># Append</span>

<span class="hljs-comment"># String processing</span>
<span class="hljs-built_in">seq</span>=<span class="hljs-string">"ATGCCGTAA"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq:0:3}</span>"</span>           <span class="hljs-comment"># Extract the first 3 characters: ATG</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq//T/U}</span>"</span>          <span class="hljs-comment"># Replace all T→U</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq/AT/XX}</span>"</span>         <span class="hljs-comment"># Replace the first AT</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq^^}</span>"</span>             <span class="hljs-comment"># Convert to uppercase (bash 4+)</span>
</code></pre>
<h2>7. Practical Scripts</h2>
<h3>7.1 Batch Rename</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># Change all .txt extensions to .log</span>
<span class="hljs-built_in">set</span> -euo pipefail

<span class="hljs-keyword">for</span> file <span class="hljs-keyword">in</span> *.txt; <span class="hljs-keyword">do</span>
    [ -e <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] || <span class="hljs-built_in">continue</span>          <span class="hljs-comment"># Skip when there is no match</span>
    <span class="hljs-built_in">mv</span> <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> <span class="hljs-string">"<span class="hljs-variable">\${file%.txt}</span>.log"</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Renamed: <span class="hljs-variable">$file</span> -> <span class="hljs-variable">\${file%.txt}</span>.log"</span>
<span class="hljs-keyword">done</span>
</code></pre>
<h3>7.2 Bioinformatics Batch Processing: Running FASTQ Quality Control</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># Usage: ./qc_pipeline.sh data_dir output_dir</span>
<span class="hljs-built_in">set</span> -euo pipefail

DATA_DIR=<span class="hljs-variable">\${1:?"Please provide the data directory"}</span>
OUT_DIR=<span class="hljs-variable">\${2:?"Please provide the output directory"}</span>
<span class="hljs-built_in">mkdir</span> -p <span class="hljs-string">"<span class="hljs-variable">$OUT_DIR</span>"</span>

<span class="hljs-keyword">for</span> fq <span class="hljs-keyword">in</span> <span class="hljs-string">"<span class="hljs-variable">$DATA_DIR</span>"</span>/*_R1.fastq.gz; <span class="hljs-keyword">do</span>
    base=$(<span class="hljs-built_in">basename</span> <span class="hljs-string">"<span class="hljs-variable">$fq</span>"</span> _R1.fastq.gz)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"=== Processing sample: <span class="hljs-variable">$base</span> ==="</span>

    fastqc -o <span class="hljs-string">"<span class="hljs-variable">$OUT_DIR</span>"</span> <span class="hljs-string">"<span class="hljs-variable">$fq</span>"</span>
    fastqc -o <span class="hljs-string">"<span class="hljs-variable">$OUT_DIR</span>"</span> <span class="hljs-string">"<span class="hljs-variable">\${DATA_DIR}</span>/<span class="hljs-variable">\${base}</span>_R2.fastq.gz"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-built_in">echo</span> <span class="hljs-string">"All done, results are in <span class="hljs-variable">$OUT_DIR</span>"</span>
</code></pre>
<h3>7.3 Log Monitoring and Alerting</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># Monitor the number of ERRORs in the log, print an alert when exceeding the threshold</span>
LOG_FILE=<span class="hljs-string">"app.log"</span>
THRESHOLD=<span class="hljs-variable">\${1:-10}</span>

count=$(grep -c <span class="hljs-string">"ERROR"</span> <span class="hljs-string">"<span class="hljs-variable">$LOG_FILE</span>"</span> || <span class="hljs-literal">true</span>)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Current ERROR count: <span class="hljs-variable">$count</span>"</span>

<span class="hljs-keyword">if</span> [ <span class="hljs-string">"<span class="hljs-variable">$count</span>"</span> -gt <span class="hljs-string">"<span class="hljs-variable">$THRESHOLD</span>"</span> ]; <span class="hljs-keyword">then</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"⚠️ Error count exceeds the threshold <span class="hljs-variable">$THRESHOLD</span>!"</span>
    <span class="hljs-built_in">exit</span> 1
<span class="hljs-keyword">fi</span>
</code></pre>
<h2>8. Common Pitfalls to Watch Out For</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. No spaces around the = sign in variable assignment</span>
name = <span class="hljs-string">"x"</span>    <span class="hljs-comment"># Wrong! Will be treated as a command</span>
name=<span class="hljs-string">"x"</span>      <span class="hljs-comment"># Correct</span>

<span class="hljs-comment"># 2. Quoting: paths containing spaces must be quoted</span>
file=<span class="hljs-string">"my data.txt"</span>
<span class="hljs-built_in">cat</span> <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span>   <span class="hljs-comment"># Correct</span>
<span class="hljs-built_in">cat</span> <span class="hljs-variable">$file</span>     <span class="hljs-comment"># Wrong: splits into two arguments</span>

<span class="hljs-comment"># 3. There must be spaces on both sides inside [ ]</span>
[ <span class="hljs-string">"<span class="hljs-variable">$a</span>"</span> = <span class="hljs-string">"<span class="hljs-variable">$b</span>"</span> ]   <span class="hljs-comment"># Correct</span>
[<span class="hljs-string">"<span class="hljs-variable">$a</span>"</span>=<span class="hljs-string">"<span class="hljs-variable">$b</span>"</span>]       <span class="hljs-comment"># Wrong</span>

<span class="hljs-comment"># 4. The pitfall of set -e in pipelines: grep returns non-zero when no match is found, which will interrupt the script</span>
grep <span class="hljs-string">"x"</span> file || <span class="hljs-literal">true</span>    <span class="hljs-comment"># Explicitly tolerate</span>
</code></pre>
<h2>9. Summary</h2>
<ul>
<li>Script header: <code>#!/bin/bash</code> + <code>set -euo pipefail</code></li>
<li>Variables: <code>$name</code>, <code>\${#name}</code>, <code>$(cmd)</code>, <code>$((...))</code></li>
<li>Conditionals: <code>[ condition ]</code>, <code>(( arithmetic ))</code>, <code>case</code></li>
<li>Loops: <code>for file in *.txt</code>, <code>while read line</code></li>
<li>Functions: <code>local</code> variables, capturing return values with echo</li>
</ul>
<p>At this point, the programming track tutorials are complete: Python ×3, R ×3, Linux ×1, Bash ×1, from zero to hands-on practice.</p>`,Eb=`<h1>Bash 编程：变量、条件、循环与实用脚本</h1>
<p>Bash 脚本把命令行操作组织成可复用的程序，是自动化批处理的核心工具。本文从脚本基础到实战，覆盖变量、条件、循环、函数与常见任务。</p>
<h2>1. 脚本基础</h2>
<h3>1.1 脚本结构</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># 第一行 shebang：指定解释器</span>

<span class="hljs-comment"># 给脚本执行权限</span>
<span class="hljs-built_in">chmod</span> +x myscript.sh
<span class="hljs-comment"># 运行</span>
./myscript.sh
<span class="hljs-comment"># 或用 bash 直接运行（无需执行权限）</span>
bash myscript.sh
</code></pre>
<h3>1.2 调试</h3>
<pre><code class="hljs language-bash">bash -x script.sh    <span class="hljs-comment"># 跟踪执行，显示每条命令</span>
<span class="hljs-comment"># 脚本内开启：set -x（跟踪）/ set -e（出错即退出）/ set -u（未定义变量报错）</span>
</code></pre>
<p>生产脚本推荐在开头加：</p>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-built_in">set</span> -euo pipefail   <span class="hljs-comment"># 出错即停、未定义变量报错、管道失败检测</span>
</code></pre>
<h2>2. 变量</h2>
<h3>2.1 定义与引用</h3>
<pre><code class="hljs language-bash">name=<span class="hljs-string">"zorrooz"</span>        <span class="hljs-comment"># 注意：= 两边不能有空格</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Hello, <span class="hljs-variable">$name</span>"</span>   <span class="hljs-comment"># 双引号内变量展开</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">'Hello, $name'</span>   <span class="hljs-comment"># 单引号内原样输出</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"长度: <span class="hljs-variable">\${#name}</span>"</span> <span class="hljs-comment"># 字符串长度</span>
</code></pre>
<h3>2.2 命令替换与算术</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 命令替换：$(...)</span>
today=$(<span class="hljs-built_in">date</span> +%F)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"今天是 <span class="hljs-variable">$today</span>"</span>
files=$(<span class="hljs-built_in">ls</span> | <span class="hljs-built_in">wc</span> -l)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"共 <span class="hljs-variable">$files</span> 个文件"</span>

<span class="hljs-comment"># 算术：$((...))</span>
a=10
b=3
<span class="hljs-built_in">echo</span> $((a + b))    <span class="hljs-comment"># 13</span>
<span class="hljs-built_in">echo</span> $((a % b))    <span class="hljs-comment"># 1</span>
</code></pre>
<h3>2.3 位置参数</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"脚本名: <span class="hljs-variable">$0</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"第一个参数: <span class="hljs-variable">$1</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"第二个参数: <span class="hljs-variable">$2</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"所有参数: <span class="hljs-variable">$@</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"参数个数: <span class="hljs-variable">$#</span>"</span>

<span class="hljs-comment"># 运行：./script.sh file1 file2</span>
</code></pre>
<h3>2.4 特殊变量</h3>
<pre><code class="hljs language-bash">$?    <span class="hljs-comment"># 上一条命令的退出码（0 成功，非 0 失败）</span>
$$    <span class="hljs-comment"># 当前进程 PID</span>
$!    <span class="hljs-comment"># 后台进程 PID</span>
</code></pre>
<h2>3. 条件判断</h2>
<h3>3.1 if / elif / else</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
score=<span class="hljs-variable">$1</span>

<span class="hljs-keyword">if</span> (( score >= <span class="hljs-number">90</span> )); <span class="hljs-keyword">then</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"优秀"</span>
<span class="hljs-keyword">elif</span> (( score >= <span class="hljs-number">60</span> )); <span class="hljs-keyword">then</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"及格"</span>
<span class="hljs-keyword">else</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"不及格"</span>
<span class="hljs-keyword">fi</span>
</code></pre>
<h3>3.2 test 表达式</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 文件判断</span>
[ -f file.txt ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"是普通文件"</span>
[ -d <span class="hljs-built_in">dir</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"是目录"</span>
[ -e path ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"存在"</span>
[ -s file ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"非空文件"</span>

<span class="hljs-comment"># 字符串比较</span>
[ <span class="hljs-string">"<span class="hljs-variable">$name</span>"</span> = <span class="hljs-string">"zorrooz"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"匹配"</span>
[ -z <span class="hljs-string">"<span class="hljs-variable">$var</span>"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"变量为空"</span>     <span class="hljs-comment"># -n 非空</span>

<span class="hljs-comment"># 数值比较（[ ] 内）</span>
[ <span class="hljs-string">"<span class="hljs-variable">$a</span>"</span> -gt <span class="hljs-string">"<span class="hljs-variable">$b</span>"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"a > b"</span>    <span class="hljs-comment"># -eq -ne -gt -ge -lt -le</span>

<span class="hljs-comment"># 逻辑组合</span>
[ -f <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] &#x26;&#x26; [ -s <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"非空文件"</span>
[ -f <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] || <span class="hljs-built_in">echo</span> <span class="hljs-string">"文件不存在"</span>
</code></pre>
<h3>3.3 case 分支</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-keyword">case</span> <span class="hljs-string">"<span class="hljs-variable">$1</span>"</span> <span class="hljs-keyword">in</span>
  start)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"启动服务"</span>
    ;;
  stop|<span class="hljs-built_in">kill</span>)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"停止服务"</span>
    ;;
  *)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"用法: <span class="hljs-variable">$0</span> {start|stop}"</span>
    <span class="hljs-built_in">exit</span> 1
    ;;
<span class="hljs-keyword">esac</span>
</code></pre>
<h2>4. 循环</h2>
<h3>4.1 for 循环</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 遍历列表</span>
<span class="hljs-keyword">for</span> fruit <span class="hljs-keyword">in</span> apple banana cherry; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"处理 <span class="hljs-variable">$fruit</span>"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-comment"># 遍历范围</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> {1..5}; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"第 <span class="hljs-variable">$i</span> 次"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-comment"># 遍历文件（最常用）</span>
<span class="hljs-keyword">for</span> file <span class="hljs-keyword">in</span> *.fasta; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"处理 <span class="hljs-variable">$file</span>"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-comment"># C 风格</span>
<span class="hljs-keyword">for</span> ((i = <span class="hljs-number">0</span>; i &#x3C; <span class="hljs-number">10</span>; i++)); <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-variable">$i</span>
<span class="hljs-keyword">done</span>
</code></pre>
<h3>4.2 while 循环</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 逐行读取文件</span>
<span class="hljs-keyword">while</span> IFS= <span class="hljs-built_in">read</span> -r line; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"读取: <span class="hljs-variable">$line</span>"</span>
<span class="hljs-keyword">done</span> &#x3C; input.txt

<span class="hljs-comment"># 计数循环</span>
count=0
<span class="hljs-keyword">while</span> [ <span class="hljs-variable">$count</span> -lt 5 ]; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-variable">$count</span>
    count=$((count + <span class="hljs-number">1</span>))
<span class="hljs-keyword">done</span>
</code></pre>
<h3>4.3 break / continue</h3>
<pre><code class="hljs language-bash"><span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> {1..10}; <span class="hljs-keyword">do</span>
    [ <span class="hljs-variable">$i</span> -eq 3 ] &#x26;&#x26; <span class="hljs-built_in">continue</span>   <span class="hljs-comment"># 跳过 3</span>
    [ <span class="hljs-variable">$i</span> -eq 7 ] &#x26;&#x26; <span class="hljs-built_in">break</span>      <span class="hljs-comment"># 到 7 终止</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-variable">$i</span>
<span class="hljs-keyword">done</span>
</code></pre>
<h2>5. 函数</h2>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>

<span class="hljs-comment"># 定义函数</span>
<span class="hljs-function"><span class="hljs-title">greet</span></span>() {
    <span class="hljs-built_in">local</span> name=<span class="hljs-variable">$1</span>           <span class="hljs-comment"># local：局部变量</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Hello, <span class="hljs-variable">$name</span>!"</span>
    <span class="hljs-built_in">return</span> 0                <span class="hljs-comment"># 返回值（退出码）</span>
}

<span class="hljs-comment"># 调用</span>
greet <span class="hljs-string">"zorrooz"</span>

<span class="hljs-comment"># 函数返回数值（注意：不是 return，是 echo 捕获）</span>
<span class="hljs-function"><span class="hljs-title">add</span></span>() {
    <span class="hljs-built_in">echo</span> $(( <span class="hljs-variable">$1</span> + <span class="hljs-variable">$2</span> ))
}
result=$(add 3 4)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"3 + 4 = <span class="hljs-variable">$result</span>"</span>
</code></pre>
<h2>6. 数组与字符串</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 数组</span>
genes=(<span class="hljs-string">"TP53"</span> <span class="hljs-string">"BRCA1"</span> <span class="hljs-string">"EGFR"</span>)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${genes[0]}</span>"</span>          <span class="hljs-comment"># TP53</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${genes[@]}</span>"</span>          <span class="hljs-comment"># 所有元素</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${#genes[@]}</span>"</span>         <span class="hljs-comment"># 长度</span>
genes+=(<span class="hljs-string">"MYC"</span>)              <span class="hljs-comment"># 追加</span>

<span class="hljs-comment"># 字符串处理</span>
<span class="hljs-built_in">seq</span>=<span class="hljs-string">"ATGCCGTAA"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq:0:3}</span>"</span>           <span class="hljs-comment"># 截取前 3 个字符：ATG</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq//T/U}</span>"</span>          <span class="hljs-comment"># 全部替换 T→U</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq/AT/XX}</span>"</span>         <span class="hljs-comment"># 替换第一个 AT</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq^^}</span>"</span>             <span class="hljs-comment"># 转大写（bash 4+）</span>
</code></pre>
<h2>7. 实战脚本</h2>
<h3>7.1 批量重命名</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># 把所有 .txt 后缀改为 .log</span>
<span class="hljs-built_in">set</span> -euo pipefail

<span class="hljs-keyword">for</span> file <span class="hljs-keyword">in</span> *.txt; <span class="hljs-keyword">do</span>
    [ -e <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] || <span class="hljs-built_in">continue</span>          <span class="hljs-comment"># 无匹配时跳过</span>
    <span class="hljs-built_in">mv</span> <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> <span class="hljs-string">"<span class="hljs-variable">\${file%.txt}</span>.log"</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"重命名: <span class="hljs-variable">$file</span> -> <span class="hljs-variable">\${file%.txt}</span>.log"</span>
<span class="hljs-keyword">done</span>
</code></pre>
<h3>7.2 生信批处理：批量跑 FASTQ 质控</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># 用法: ./qc_pipeline.sh data_dir output_dir</span>
<span class="hljs-built_in">set</span> -euo pipefail

DATA_DIR=<span class="hljs-variable">\${1:?"请提供数据目录"}</span>
OUT_DIR=<span class="hljs-variable">\${2:?"请提供输出目录"}</span>
<span class="hljs-built_in">mkdir</span> -p <span class="hljs-string">"<span class="hljs-variable">$OUT_DIR</span>"</span>

<span class="hljs-keyword">for</span> fq <span class="hljs-keyword">in</span> <span class="hljs-string">"<span class="hljs-variable">$DATA_DIR</span>"</span>/*_R1.fastq.gz; <span class="hljs-keyword">do</span>
    base=$(<span class="hljs-built_in">basename</span> <span class="hljs-string">"<span class="hljs-variable">$fq</span>"</span> _R1.fastq.gz)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"=== 处理样本: <span class="hljs-variable">$base</span> ==="</span>

    fastqc -o <span class="hljs-string">"<span class="hljs-variable">$OUT_DIR</span>"</span> <span class="hljs-string">"<span class="hljs-variable">$fq</span>"</span>
    fastqc -o <span class="hljs-string">"<span class="hljs-variable">$OUT_DIR</span>"</span> <span class="hljs-string">"<span class="hljs-variable">\${DATA_DIR}</span>/<span class="hljs-variable">\${base}</span>_R2.fastq.gz"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-built_in">echo</span> <span class="hljs-string">"全部完成，结果位于 <span class="hljs-variable">$OUT_DIR</span>"</span>
</code></pre>
<h3>7.3 日志监控与告警</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># 监控日志中 ERROR 数量，超过阈值打印告警</span>
LOG_FILE=<span class="hljs-string">"app.log"</span>
THRESHOLD=<span class="hljs-variable">\${1:-10}</span>

count=$(grep -c <span class="hljs-string">"ERROR"</span> <span class="hljs-string">"<span class="hljs-variable">$LOG_FILE</span>"</span> || <span class="hljs-literal">true</span>)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"当前 ERROR 数量: <span class="hljs-variable">$count</span>"</span>

<span class="hljs-keyword">if</span> [ <span class="hljs-string">"<span class="hljs-variable">$count</span>"</span> -gt <span class="hljs-string">"<span class="hljs-variable">$THRESHOLD</span>"</span> ]; <span class="hljs-keyword">then</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"⚠️ 错误数量超过阈值 <span class="hljs-variable">$THRESHOLD</span>！"</span>
    <span class="hljs-built_in">exit</span> 1
<span class="hljs-keyword">fi</span>
</code></pre>
<h2>8. 常见坑位提醒</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. 变量赋值等号两边不能有空格</span>
name = <span class="hljs-string">"x"</span>    <span class="hljs-comment"># 错误！会被当成命令</span>
name=<span class="hljs-string">"x"</span>      <span class="hljs-comment"># 正确</span>

<span class="hljs-comment"># 2. 引号：路径含空格必须加引号</span>
file=<span class="hljs-string">"my data.txt"</span>
<span class="hljs-built_in">cat</span> <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span>   <span class="hljs-comment"># 正确</span>
<span class="hljs-built_in">cat</span> <span class="hljs-variable">$file</span>     <span class="hljs-comment"># 错误：拆成两个参数</span>

<span class="hljs-comment"># 3. [ ] 内部两端必须有空格</span>
[ <span class="hljs-string">"<span class="hljs-variable">$a</span>"</span> = <span class="hljs-string">"<span class="hljs-variable">$b</span>"</span> ]   <span class="hljs-comment"># 正确</span>
[<span class="hljs-string">"<span class="hljs-variable">$a</span>"</span>=<span class="hljs-string">"<span class="hljs-variable">$b</span>"</span>]       <span class="hljs-comment"># 错误</span>

<span class="hljs-comment"># 4. 管道中 set -e 的坑：grep 无匹配返回非 0 会中断</span>
grep <span class="hljs-string">"x"</span> file || <span class="hljs-literal">true</span>    <span class="hljs-comment"># 显式容忍</span>
</code></pre>
<h2>9. 小结</h2>
<ul>
<li>脚本头：<code>#!/bin/bash</code> + <code>set -euo pipefail</code></li>
<li>变量：<code>$name</code>、<code>\${#name}</code>、<code>$(cmd)</code>、<code>$((...))</code></li>
<li>判断：<code>[ 条件 ]</code>、<code>(( 算术 ))</code>、<code>case</code></li>
<li>循环：<code>for file in *.txt</code>、<code>while read line</code></li>
<li>函数：<code>local</code> 变量、echo 捕获返回值</li>
</ul>
<p>至此编程方向教程完成：Python ×3、R ×3、Linux ×1、Bash ×1，从零到实战。</p>`,Sb=`<h1>Linux Command Line Basics: File System, Permissions, and Text Processing</h1>
<p>The Linux command line (Shell) is a fundamental skill for bioinformatics and scientific computing: server operations, workflow construction, and batch processing all depend on it. This article covers the most commonly used commands, and it is recommended to practice while learning using WSL or a remote server.</p>
<h2>1. Terminal and Shell</h2>
<pre><code class="hljs language-bash"><span class="hljs-built_in">echo</span> <span class="hljs-variable">$SHELL</span>        <span class="hljs-comment"># View the current shell, usually /bin/bash</span>
<span class="hljs-built_in">whoami</span>             <span class="hljs-comment"># Current user</span>
<span class="hljs-built_in">pwd</span>                <span class="hljs-comment"># Current directory (print working directory)</span>
</code></pre>
<h2>2. File System Navigation</h2>
<pre><code class="hljs language-bash"><span class="hljs-built_in">ls</span>                  <span class="hljs-comment"># List files</span>
<span class="hljs-built_in">ls</span> -l               <span class="hljs-comment"># Detailed information (permissions, size, time)</span>
<span class="hljs-built_in">ls</span> -a               <span class="hljs-comment"># Include hidden files (starting with .)</span>
<span class="hljs-built_in">ls</span> -lh              <span class="hljs-comment"># Human-readable sizes</span>

<span class="hljs-built_in">cd</span> /home/user       <span class="hljs-comment"># Enter directory</span>
<span class="hljs-built_in">cd</span> ..               <span class="hljs-comment"># Parent directory</span>
<span class="hljs-built_in">cd</span> ~                <span class="hljs-comment"># User home directory</span>
<span class="hljs-built_in">cd</span> -                <span class="hljs-comment"># Previous directory</span>

<span class="hljs-built_in">pwd</span>                 <span class="hljs-comment"># Display current location</span>
</code></pre>
<p><strong>Path rules</strong>: <code>/</code> is the root directory; <code>.</code> is the current directory; <code>..</code> is the parent directory; <code>~</code> is the home directory.</p>
<h2>3. File and Directory Operations</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Create</span>
<span class="hljs-built_in">touch</span> file.txt      <span class="hljs-comment"># Create an empty file</span>
<span class="hljs-built_in">mkdir</span> data          <span class="hljs-comment"># Create a directory</span>
<span class="hljs-built_in">mkdir</span> -p a/b/c      <span class="hljs-comment"># Recursively create multi-level directories</span>

<span class="hljs-comment"># Copy / Move / Delete</span>
<span class="hljs-built_in">cp</span> file.txt copy.txt
<span class="hljs-built_in">cp</span> -r data/ data_backup/    <span class="hljs-comment"># Copy directory</span>
<span class="hljs-built_in">mv</span> file.txt newname.txt     <span class="hljs-comment"># Rename/move</span>
<span class="hljs-built_in">rm</span> file.txt                 <span class="hljs-comment"># Delete file</span>
<span class="hljs-built_in">rm</span> -r data/                 <span class="hljs-comment"># Delete directory (dangerous!)</span>
<span class="hljs-built_in">rm</span> -rf data/                <span class="hljs-comment"># Force recursive delete (use with caution)</span>

<span class="hljs-comment"># View contents</span>
<span class="hljs-built_in">cat</span> file.txt        <span class="hljs-comment"># Output all contents</span>
less file.txt       <span class="hljs-comment"># Page through (q to quit, / to search)</span>
<span class="hljs-built_in">head</span> -n 5 file.txt  <span class="hljs-comment"># First 5 lines</span>
<span class="hljs-built_in">tail</span> -n 5 file.txt  <span class="hljs-comment"># Last 5 lines</span>
<span class="hljs-built_in">tail</span> -f log.txt     <span class="hljs-comment"># Follow in real time (commonly used for viewing logs)</span>
</code></pre>
<blockquote>
<p><strong><code>rm -rf</code> is a dangerous command</strong>: confirm the path before executing, or first <code>ls</code> to verify.</p>
</blockquote>
<h2>4. Wildcards and Redirection</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Wildcards</span>
<span class="hljs-built_in">ls</span> *.fasta          <span class="hljs-comment"># All .fasta files</span>
<span class="hljs-built_in">ls</span> data_??.txt      <span class="hljs-comment"># ? matches a single character</span>
<span class="hljs-built_in">ls</span> [abc]*           <span class="hljs-comment"># Starts with a/b/c</span>

<span class="hljs-comment"># Redirection</span>
<span class="hljs-built_in">ls</span> > list.txt       <span class="hljs-comment"># Overwrite</span>
<span class="hljs-built_in">ls</span> >> list.txt      <span class="hljs-comment"># Append</span>
<span class="hljs-built_in">ls</span> 2> error.log     <span class="hljs-comment"># Redirect error output</span>
<span class="hljs-built_in">ls</span> 2>/dev/null      <span class="hljs-comment"># Discard error output</span>

<span class="hljs-comment"># Pipe: output of the previous command becomes input of the next</span>
<span class="hljs-built_in">ls</span> -l | <span class="hljs-built_in">wc</span> -l                   <span class="hljs-comment"># Count number of files</span>
<span class="hljs-built_in">history</span> | grep python           <span class="hljs-comment"># Search history commands</span>
</code></pre>
<h2>5. Permission Management</h2>
<pre><code class="hljs language-bash"><span class="hljs-built_in">ls</span> -l
<span class="hljs-comment"># -rw-r--r-- 1 user group 1234 Aug 4 10:00 file.txt</span>
<span class="hljs-comment"># ^permissions       ^owner ^group</span>

<span class="hljs-built_in">chmod</span> 755 script.sh    <span class="hljs-comment"># rwxr-xr-x: owner read/write/execute, others read/execute only</span>
<span class="hljs-built_in">chmod</span> +x script.sh     <span class="hljs-comment"># Add execute permission (required to run scripts)</span>
<span class="hljs-built_in">chmod</span> -w file.txt      <span class="hljs-comment"># Remove write permission</span>

<span class="hljs-built_in">chown</span> user:group file  <span class="hljs-comment"># Change owner and group (requires root)</span>
</code></pre>
<p>Permission numbers: <code>r=4</code>, <code>w=2</code>, <code>x=1</code>, and the three digits represent owner/group/others respectively.</p>
<h2>6. The Three Musketeers of Text Processing</h2>
<h3>6.1 grep: Search</h3>
<pre><code class="hljs language-bash">grep <span class="hljs-string">"TP53"</span> genes.txt          <span class="hljs-comment"># Find lines containing TP53</span>
grep -i <span class="hljs-string">"tp53"</span> genes.txt       <span class="hljs-comment"># Ignore case</span>
grep -v <span class="hljs-string">"comment"</span> file.txt     <span class="hljs-comment"># Reverse match (exclude)</span>
grep -c <span class="hljs-string">"gene"</span> file.txt        <span class="hljs-comment"># Count</span>
grep -n <span class="hljs-string">"pattern"</span> file.txt     <span class="hljs-comment"># Show line numbers</span>
grep -r <span class="hljs-string">"TODO"</span> src/            <span class="hljs-comment"># Recursively search directory</span>

<span class="hljs-comment"># Pipe combinations</span>
ps aux | grep python           <span class="hljs-comment"># Find python processes</span>
<span class="hljs-built_in">cat</span> reads.fastq | grep -c <span class="hljs-string">"^@"</span> <span class="hljs-comment"># Count number of sequences in FASTQ</span>
</code></pre>
<h3>6.2 sed: Stream Editing</h3>
<pre><code class="hljs language-bash">sed -n <span class="hljs-string">'5,10p'</span> file.txt        <span class="hljs-comment"># Print lines 5-10</span>
sed <span class="hljs-string">'s/old/new/'</span> file.txt      <span class="hljs-comment"># Replace first match in each line</span>
sed <span class="hljs-string">'s/old/new/g'</span> file.txt     <span class="hljs-comment"># Global replace</span>
sed -i <span class="hljs-string">'s/old/new/g'</span> file.txt  <span class="hljs-comment"># Modify file in place (-i in-place)</span>
sed <span class="hljs-string">'/^#/d'</span> config.conf        <span class="hljs-comment"># Delete comment lines</span>
</code></pre>
<h3>6.3 awk: Column-Based Processing</h3>
<pre><code class="hljs language-bash">awk <span class="hljs-string">'{print $1}'</span> file.txt      <span class="hljs-comment"># Print first column</span>
awk -F<span class="hljs-string">','</span> <span class="hljs-string">'{print $1, $3}'</span> data.csv   <span class="hljs-comment"># Separate by comma</span>
awk <span class="hljs-string">'$3 > 50 {print $1}'</span> data.txt     <span class="hljs-comment"># Conditional filter</span>
awk <span class="hljs-string">'{sum += $2} END {print sum}'</span> data.txt   <span class="hljs-comment"># Sum second column</span>
awk <span class="hljs-string">'NR > 1 {print}'</span> file.txt  <span class="hljs-comment"># Skip header line</span>
</code></pre>
<blockquote>
<p>awk field separator defaults to whitespace; <code>-F','</code> specifies comma.</p>
</blockquote>
<h2>7. Process Management</h2>
<pre><code class="hljs language-bash">ps aux                <span class="hljs-comment"># View all processes</span>
ps aux | grep python  <span class="hljs-comment"># Find specific processes</span>
top                   <span class="hljs-comment"># Real-time monitoring (q to quit)</span>
htop                  <span class="hljs-comment"># Enhanced monitoring (more intuitive)</span>

<span class="hljs-built_in">kill</span> PID              <span class="hljs-comment"># Terminate process</span>
<span class="hljs-built_in">kill</span> -9 PID           <span class="hljs-comment"># Force kill</span>
pkill python          <span class="hljs-comment"># Terminate by name</span>

<span class="hljs-comment"># Background execution</span>
python script.py &#x26;    <span class="hljs-comment"># Run in background</span>
<span class="hljs-built_in">nohup</span> python script.py > log.txt 2>&#x26;1 &#x26;   <span class="hljs-comment"># Run detached from terminal (common on servers)</span>
<span class="hljs-built_in">jobs</span>                  <span class="hljs-comment"># View background jobs</span>
<span class="hljs-built_in">fg</span>                    <span class="hljs-comment"># Bring to foreground</span>
</code></pre>
<blockquote>
<p>For long-running tasks on remote servers, use <code>nohup ... &#x26;</code> or <code>tmux</code> / <code>screen</code> to avoid task interruption due to terminal disconnection.</p>
</blockquote>
<h2>8. Practical Combinations</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Count files and sizes</span>
<span class="hljs-built_in">ls</span> -lh | awk <span class="hljs-string">'{print $5, $9}'</span>
<span class="hljs-built_in">du</span> -sh data/          <span class="hljs-comment"># Total size of directory</span>
<span class="hljs-built_in">df</span> -h                 <span class="hljs-comment"># Disk space</span>

<span class="hljs-comment"># Compression and decompression</span>
tar -czf archive.tar.gz <span class="hljs-built_in">dir</span>/    <span class="hljs-comment"># Pack and compress</span>
tar -xzf archive.tar.gz         <span class="hljs-comment"># Decompress</span>
zip -r archive.zip <span class="hljs-built_in">dir</span>/
unzip archive.zip

<span class="hljs-comment"># Find files</span>
find . -name <span class="hljs-string">"*.log"</span>            <span class="hljs-comment"># Find by name</span>
find . -size +100M              <span class="hljs-comment"># Find by size</span>
<span class="hljs-built_in">which</span> python                    <span class="hljs-comment"># Path of the command</span>
</code></pre>
<h2>9. WSL and Development Environment</h2>
<p>Under Windows, WSL (Windows Subsystem for Linux) is recommended:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Install (in administrator PowerShell)</span>
wsl --install

<span class="hljs-comment"># Enter the Ubuntu subsystem</span>
wsl

<span class="hljs-comment"># Configure the development environment in WSL</span>
<span class="hljs-built_in">sudo</span> apt update
<span class="hljs-built_in">sudo</span> apt install build-essential git python3 python3-pip
</code></pre>
<p>VS Code connects to WSL: after installing the Remote - WSL extension, <code>code .</code> opens directly.</p>
<h2>10. Summary</h2>
<ul>
<li>Navigation: <code>cd</code> / <code>ls</code> / <code>pwd</code>; operations: <code>cp</code> / <code>mv</code> / <code>rm</code> / <code>mkdir</code></li>
<li>The three musketeers: <code>grep</code> for searching, <code>sed</code> for replacing, <code>awk</code> for column-based processing</li>
<li>The pipe <code>|</code> and redirection <code>></code> are the core of combining commands</li>
<li>Permissions: <code>chmod 755</code> / <code>chmod +x</code>; processes: <code>ps</code> / <code>kill</code> / <code>nohup</code></li>
</ul>
<p>The next article will introduce Bash programming: organizing commands into reusable scripts.</p>`,Tb=`<h1>Linux 命令行基础：文件系统、权限与文本处理</h1>
<p>Linux 命令行（Shell）是生物信息学与科学计算的基本功：服务器操作、流程搭建、批量处理都离不开它。本文覆盖最常用的命令，建议配合 WSL 或远程服务器边学边练。</p>
<h2>1. 终端与 Shell</h2>
<pre><code class="hljs language-bash"><span class="hljs-built_in">echo</span> <span class="hljs-variable">$SHELL</span>        <span class="hljs-comment"># 查看当前 shell，通常为 /bin/bash</span>
<span class="hljs-built_in">whoami</span>             <span class="hljs-comment"># 当前用户</span>
<span class="hljs-built_in">pwd</span>                <span class="hljs-comment"># 当前目录（print working directory）</span>
</code></pre>
<h2>2. 文件系统导航</h2>
<pre><code class="hljs language-bash"><span class="hljs-built_in">ls</span>                  <span class="hljs-comment"># 列出文件</span>
<span class="hljs-built_in">ls</span> -l               <span class="hljs-comment"># 详细信息（权限、大小、时间）</span>
<span class="hljs-built_in">ls</span> -a               <span class="hljs-comment"># 包含隐藏文件（以 . 开头）</span>
<span class="hljs-built_in">ls</span> -lh              <span class="hljs-comment"># 人类可读大小</span>

<span class="hljs-built_in">cd</span> /home/user       <span class="hljs-comment"># 进入目录</span>
<span class="hljs-built_in">cd</span> ..               <span class="hljs-comment"># 上级目录</span>
<span class="hljs-built_in">cd</span> ~                <span class="hljs-comment"># 用户主目录</span>
<span class="hljs-built_in">cd</span> -                <span class="hljs-comment"># 上一个目录</span>

<span class="hljs-built_in">pwd</span>                 <span class="hljs-comment"># 显示当前位置</span>
</code></pre>
<p><strong>路径规则</strong>：<code>/</code> 是根目录；<code>.</code> 当前目录；<code>..</code> 上级目录；<code>~</code> 主目录。</p>
<h2>3. 文件与目录操作</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 创建</span>
<span class="hljs-built_in">touch</span> file.txt      <span class="hljs-comment"># 创建空文件</span>
<span class="hljs-built_in">mkdir</span> data          <span class="hljs-comment"># 创建目录</span>
<span class="hljs-built_in">mkdir</span> -p a/b/c      <span class="hljs-comment"># 递归创建多级目录</span>

<span class="hljs-comment"># 复制 / 移动 / 删除</span>
<span class="hljs-built_in">cp</span> file.txt copy.txt
<span class="hljs-built_in">cp</span> -r data/ data_backup/    <span class="hljs-comment"># 复制目录</span>
<span class="hljs-built_in">mv</span> file.txt newname.txt     <span class="hljs-comment"># 重命名/移动</span>
<span class="hljs-built_in">rm</span> file.txt                 <span class="hljs-comment"># 删除文件</span>
<span class="hljs-built_in">rm</span> -r data/                 <span class="hljs-comment"># 删除目录（危险！）</span>
<span class="hljs-built_in">rm</span> -rf data/                <span class="hljs-comment"># 强制递归删除（慎用）</span>

<span class="hljs-comment"># 查看内容</span>
<span class="hljs-built_in">cat</span> file.txt        <span class="hljs-comment"># 输出全部</span>
less file.txt       <span class="hljs-comment"># 分页浏览（q 退出，/搜索）</span>
<span class="hljs-built_in">head</span> -n 5 file.txt  <span class="hljs-comment"># 前 5 行</span>
<span class="hljs-built_in">tail</span> -n 5 file.txt  <span class="hljs-comment"># 后 5 行</span>
<span class="hljs-built_in">tail</span> -f log.txt     <span class="hljs-comment"># 实时跟踪（查看日志常用）</span>
</code></pre>
<blockquote>
<p><strong>rm -rf 是危险命令</strong>：执行前确认路径，或先 <code>ls</code> 验证。</p>
</blockquote>
<h2>4. 通配符与重定向</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 通配符</span>
<span class="hljs-built_in">ls</span> *.fasta          <span class="hljs-comment"># 所有 .fasta 文件</span>
<span class="hljs-built_in">ls</span> data_??.txt      <span class="hljs-comment"># ? 匹配单个字符</span>
<span class="hljs-built_in">ls</span> [abc]*           <span class="hljs-comment"># 以 a/b/c 开头</span>

<span class="hljs-comment"># 重定向</span>
<span class="hljs-built_in">ls</span> > list.txt       <span class="hljs-comment"># 覆盖写入</span>
<span class="hljs-built_in">ls</span> >> list.txt      <span class="hljs-comment"># 追加写入</span>
<span class="hljs-built_in">ls</span> 2> error.log     <span class="hljs-comment"># 错误输出重定向</span>
<span class="hljs-built_in">ls</span> 2>/dev/null      <span class="hljs-comment"># 丢弃错误输出</span>

<span class="hljs-comment"># 管道：前一个命令的输出作为后一个的输入</span>
<span class="hljs-built_in">ls</span> -l | <span class="hljs-built_in">wc</span> -l                   <span class="hljs-comment"># 统计文件数</span>
<span class="hljs-built_in">history</span> | grep python           <span class="hljs-comment"># 在历史命令中搜索</span>
</code></pre>
<h2>5. 权限管理</h2>
<pre><code class="hljs language-bash"><span class="hljs-built_in">ls</span> -l
<span class="hljs-comment"># -rw-r--r-- 1 user group 1234 Aug 4 10:00 file.txt</span>
<span class="hljs-comment"># ^权限       ^属主 ^属组</span>

<span class="hljs-built_in">chmod</span> 755 script.sh    <span class="hljs-comment"># rwxr-xr-x：属主读写执行，其他只读执行</span>
<span class="hljs-built_in">chmod</span> +x script.sh     <span class="hljs-comment"># 添加执行权限（运行脚本必需）</span>
<span class="hljs-built_in">chmod</span> -w file.txt      <span class="hljs-comment"># 去除写权限</span>

<span class="hljs-built_in">chown</span> user:group file  <span class="hljs-comment"># 修改属主属组（需要 root）</span>
</code></pre>
<p>权限数字：<code>r=4</code>、<code>w=2</code>、<code>x=1</code>，三位分别表示属主/属组/其他。</p>
<h2>6. 文本处理三剑客</h2>
<h3>6.1 grep：搜索</h3>
<pre><code class="hljs language-bash">grep <span class="hljs-string">"TP53"</span> genes.txt          <span class="hljs-comment"># 查找包含 TP53 的行</span>
grep -i <span class="hljs-string">"tp53"</span> genes.txt       <span class="hljs-comment"># 忽略大小写</span>
grep -v <span class="hljs-string">"comment"</span> file.txt     <span class="hljs-comment"># 反向匹配（排除）</span>
grep -c <span class="hljs-string">"gene"</span> file.txt        <span class="hljs-comment"># 计数</span>
grep -n <span class="hljs-string">"pattern"</span> file.txt     <span class="hljs-comment"># 显示行号</span>
grep -r <span class="hljs-string">"TODO"</span> src/            <span class="hljs-comment"># 递归搜索目录</span>

<span class="hljs-comment"># 管道组合</span>
ps aux | grep python           <span class="hljs-comment"># 找 python 进程</span>
<span class="hljs-built_in">cat</span> reads.fastq | grep -c <span class="hljs-string">"^@"</span> <span class="hljs-comment"># 统计 FASTQ 中序列条数</span>
</code></pre>
<h3>6.2 sed：流编辑</h3>
<pre><code class="hljs language-bash">sed -n <span class="hljs-string">'5,10p'</span> file.txt        <span class="hljs-comment"># 打印第 5-10 行</span>
sed <span class="hljs-string">'s/old/new/'</span> file.txt      <span class="hljs-comment"># 替换每行第一个匹配</span>
sed <span class="hljs-string">'s/old/new/g'</span> file.txt     <span class="hljs-comment"># 全局替换</span>
sed -i <span class="hljs-string">'s/old/new/g'</span> file.txt  <span class="hljs-comment"># 直接修改文件（-i 原地）</span>
sed <span class="hljs-string">'/^#/d'</span> config.conf        <span class="hljs-comment"># 删除注释行</span>
</code></pre>
<h3>6.3 awk：按列处理</h3>
<pre><code class="hljs language-bash">awk <span class="hljs-string">'{print $1}'</span> file.txt      <span class="hljs-comment"># 打印第一列</span>
awk -F<span class="hljs-string">','</span> <span class="hljs-string">'{print $1, $3}'</span> data.csv   <span class="hljs-comment"># 按逗号分隔</span>
awk <span class="hljs-string">'$3 > 50 {print $1}'</span> data.txt     <span class="hljs-comment"># 条件过滤</span>
awk <span class="hljs-string">'{sum += $2} END {print sum}'</span> data.txt   <span class="hljs-comment"># 第二列求和</span>
awk <span class="hljs-string">'NR > 1 {print}'</span> file.txt  <span class="hljs-comment"># 跳过表头</span>
</code></pre>
<blockquote>
<p>awk 字段分隔符默认空白；<code>-F','</code> 指定逗号。</p>
</blockquote>
<h2>7. 进程管理</h2>
<pre><code class="hljs language-bash">ps aux                <span class="hljs-comment"># 查看所有进程</span>
ps aux | grep python  <span class="hljs-comment"># 查找特定进程</span>
top                   <span class="hljs-comment"># 实时监控（q 退出）</span>
htop                  <span class="hljs-comment"># 增强版监控（更直观）</span>

<span class="hljs-built_in">kill</span> PID              <span class="hljs-comment"># 终止进程</span>
<span class="hljs-built_in">kill</span> -9 PID           <span class="hljs-comment"># 强制终止</span>
pkill python          <span class="hljs-comment"># 按名称终止</span>

<span class="hljs-comment"># 后台运行</span>
python script.py &#x26;    <span class="hljs-comment"># 后台运行</span>
<span class="hljs-built_in">nohup</span> python script.py > log.txt 2>&#x26;1 &#x26;   <span class="hljs-comment"># 脱离终端运行（服务器常用）</span>
<span class="hljs-built_in">jobs</span>                  <span class="hljs-comment"># 查看后台任务</span>
<span class="hljs-built_in">fg</span>                    <span class="hljs-comment"># 调回前台</span>
</code></pre>
<blockquote>
<p>远程服务器跑长任务用 <code>nohup ... &#x26;</code> 或 <code>tmux</code> / <code>screen</code>，避免终端断开导致任务中断。</p>
</blockquote>
<h2>8. 实用组合</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 统计文件数量与大小</span>
<span class="hljs-built_in">ls</span> -lh | awk <span class="hljs-string">'{print $5, $9}'</span>
<span class="hljs-built_in">du</span> -sh data/          <span class="hljs-comment"># 目录总大小</span>
<span class="hljs-built_in">df</span> -h                 <span class="hljs-comment"># 磁盘空间</span>

<span class="hljs-comment"># 压缩与解压</span>
tar -czf archive.tar.gz <span class="hljs-built_in">dir</span>/    <span class="hljs-comment"># 打包压缩</span>
tar -xzf archive.tar.gz         <span class="hljs-comment"># 解压</span>
zip -r archive.zip <span class="hljs-built_in">dir</span>/
unzip archive.zip

<span class="hljs-comment"># 查找文件</span>
find . -name <span class="hljs-string">"*.log"</span>            <span class="hljs-comment"># 按名查找</span>
find . -size +100M              <span class="hljs-comment"># 按大小查找</span>
<span class="hljs-built_in">which</span> python                    <span class="hljs-comment"># 命令所在路径</span>
</code></pre>
<h2>9. WSL 与开发环境</h2>
<p>Windows 下推荐 WSL（Windows Subsystem for Linux）：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 安装（管理员 PowerShell）</span>
wsl --install

<span class="hljs-comment"># 进入 Ubuntu 子系统</span>
wsl

<span class="hljs-comment"># 在 WSL 中配置开发环境</span>
<span class="hljs-built_in">sudo</span> apt update
<span class="hljs-built_in">sudo</span> apt install build-essential git python3 python3-pip
</code></pre>
<p>VS Code 连接 WSL：安装 Remote - WSL 扩展后，<code>code .</code> 直接打开。</p>
<h2>10. 小结</h2>
<ul>
<li>导航：<code>cd</code> / <code>ls</code> / <code>pwd</code>；操作：<code>cp</code> / <code>mv</code> / <code>rm</code> / <code>mkdir</code></li>
<li>三剑客：<code>grep</code> 搜索、<code>sed</code> 替换、<code>awk</code> 按列处理</li>
<li>管道 <code>|</code> 与重定向 <code>></code> 是组合命令的核心</li>
<li>权限 <code>chmod 755</code> / <code>chmod +x</code>；进程 <code>ps</code> / <code>kill</code> / <code>nohup</code></li>
</ul>
<p>下一篇将介绍 Bash 编程：把命令组织成可复用的脚本。</p>`,Ab=`<h1>Python Advanced: Functions, Classes, and Modules</h1>
<p>After mastering the basics, this article dives into Python's functional and object-oriented programming, helping you write clean, reusable code.</p>
<h2>1. Functions</h2>
<h3>1.1 Definition and Call</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">greet</span>(<span class="hljs-params">name, greeting=<span class="hljs-string">"Hello"</span></span>):
    <span class="hljs-string">"""Return a greeting. greeting is a parameter with a default value."""</span>
    <span class="hljs-keyword">return</span> <span class="hljs-string">f"<span class="hljs-subst">{greeting}</span>, <span class="hljs-subst">{name}</span>!"</span>

<span class="hljs-built_in">print</span>(greet(<span class="hljs-string">"zorrooz"</span>))            <span class="hljs-comment"># Hello, zorrooz!</span>
<span class="hljs-built_in">print</span>(greet(<span class="hljs-string">"zorrooz"</span>, <span class="hljs-string">"Hi"</span>))      <span class="hljs-comment"># Hi, zorrooz!</span>
</code></pre>
<h3>1.2 Parameter Passing</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">show</span>(<span class="hljs-params">a, b, *args, **kwargs</span>):
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"a=<span class="hljs-subst">{a}</span>, b=<span class="hljs-subst">{b}</span>"</span>)
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"Positional args: <span class="hljs-subst">{args}</span>"</span>)      <span class="hljs-comment"># tuple</span>
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"Keyword args: <span class="hljs-subst">{kwargs}</span>"</span>)  <span class="hljs-comment"># dict</span>

show(<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>, x=<span class="hljs-number">5</span>, y=<span class="hljs-number">6</span>)
<span class="hljs-comment"># a=1, b=2</span>
<span class="hljs-comment"># Positional args: (3, 4)</span>
<span class="hljs-comment"># Keyword args: {'x': 5, 'y': 6}</span>
</code></pre>
<ul>
<li><code>*args</code>: collects extra positional arguments into a tuple</li>
<li><code>**kwargs</code>: collects extra keyword arguments into a dictionary</li>
</ul>
<h3>1.3 Unpacking Arguments</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">add</span>(<span class="hljs-params">a, b, c</span>):
    <span class="hljs-keyword">return</span> a + b + c

nums = [<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>]
<span class="hljs-built_in">print</span>(add(*nums))            <span class="hljs-comment"># 6, list unpacking</span>

data = {<span class="hljs-string">"a"</span>: <span class="hljs-number">1</span>, <span class="hljs-string">"b"</span>: <span class="hljs-number">2</span>, <span class="hljs-string">"c"</span>: <span class="hljs-number">3</span>}
<span class="hljs-built_in">print</span>(add(**data))           <span class="hljs-comment"># 6, dict unpacking</span>
</code></pre>
<h3>1.4 Lambda Anonymous Functions</h3>
<pre><code class="hljs language-python">square = <span class="hljs-keyword">lambda</span> x: x ** <span class="hljs-number">2</span>
<span class="hljs-built_in">print</span>(square(<span class="hljs-number">5</span>))             <span class="hljs-comment"># 25</span>

<span class="hljs-comment"># Combined with sorted / map / filter</span>
words = [<span class="hljs-string">"banana"</span>, <span class="hljs-string">"apple"</span>, <span class="hljs-string">"cherry"</span>]
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">sorted</span>(words, key=<span class="hljs-keyword">lambda</span> w: <span class="hljs-built_in">len</span>(w)))
<span class="hljs-comment"># ['apple', 'banana', 'cherry']</span>

nums = [<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>, <span class="hljs-number">5</span>]
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">list</span>(<span class="hljs-built_in">map</span>(<span class="hljs-keyword">lambda</span> x: x * <span class="hljs-number">2</span>, nums)))       <span class="hljs-comment"># [2, 4, 6, 8, 10]</span>
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">list</span>(<span class="hljs-built_in">filter</span>(<span class="hljs-keyword">lambda</span> x: x % <span class="hljs-number">2</span> == <span class="hljs-number">0</span>, nums)))  <span class="hljs-comment"># [2, 4]</span>
</code></pre>
<h3>1.5 Closures</h3>
<p>Define a function inside another function and reference outer variables:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">make_counter</span>():
    count = <span class="hljs-number">0</span>
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">counter</span>():
        <span class="hljs-keyword">nonlocal</span> count
        count += <span class="hljs-number">1</span>
        <span class="hljs-keyword">return</span> count
    <span class="hljs-keyword">return</span> counter

c = make_counter()
<span class="hljs-built_in">print</span>(c())   <span class="hljs-comment"># 1</span>
<span class="hljs-built_in">print</span>(c())   <span class="hljs-comment"># 2</span>
</code></pre>
<p><code>nonlocal</code> is used to modify variables in the outer function.</p>
<h2>2. Decorators</h2>
<p>A decorator is a higher-order function that "takes a function and returns a function," used to enhance its behavior without modifying the original function:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> time

<span class="hljs-keyword">def</span> <span class="hljs-title function_">timer</span>(<span class="hljs-params">func</span>):
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">wrapper</span>(<span class="hljs-params">*args, **kwargs</span>):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{func.__name__}</span> took <span class="hljs-subst">{elapsed:<span class="hljs-number">.4</span>f}</span>s"</span>)
        <span class="hljs-keyword">return</span> result
    <span class="hljs-keyword">return</span> wrapper

<span class="hljs-meta">@timer</span>
<span class="hljs-keyword">def</span> <span class="hljs-title function_">slow_sum</span>(<span class="hljs-params">n</span>):
    <span class="hljs-keyword">return</span> <span class="hljs-built_in">sum</span>(<span class="hljs-built_in">range</span>(n))

<span class="hljs-built_in">print</span>(slow_sum(<span class="hljs-number">10_000_000</span>))
<span class="hljs-comment"># slow_sum took 0.35xx s</span>
</code></pre>
<p>Decorator with parameters:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">repeat</span>(<span class="hljs-params">times</span>):
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">decorator</span>(<span class="hljs-params">func</span>):
        <span class="hljs-keyword">def</span> <span class="hljs-title function_">wrapper</span>(<span class="hljs-params">*args, **kwargs</span>):
            <span class="hljs-keyword">for</span> _ <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(times):
                func(*args, **kwargs)
        <span class="hljs-keyword">return</span> wrapper
    <span class="hljs-keyword">return</span> decorator

<span class="hljs-meta">@repeat(<span class="hljs-params"><span class="hljs-number">3</span></span>)</span>
<span class="hljs-keyword">def</span> <span class="hljs-title function_">hello</span>():
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"Hi!"</span>)

hello()
<span class="hljs-comment"># Hi! Hi! Hi!</span>
</code></pre>
<h2>3. Classes and Object-Oriented Programming</h2>
<h3>3.1 Basic Definition</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Sequence</span>:
    <span class="hljs-string">"""Represents a biological sequence."""</span>

    <span class="hljs-comment"># Class attribute: shared by all instances</span>
    alphabet = <span class="hljs-string">"ACGT"</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, seq_id, seq</span>):
        <span class="hljs-string">"""Constructor: initialize instance attributes."""</span>
        <span class="hljs-variable language_">self</span>.seq_id = seq_id
        <span class="hljs-variable language_">self</span>.seq = seq.upper()

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">length</span>(<span class="hljs-params">self</span>):
        <span class="hljs-keyword">return</span> <span class="hljs-built_in">len</span>(<span class="hljs-variable language_">self</span>.seq)

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">gc_content</span>(<span class="hljs-params">self</span>):
        gc = <span class="hljs-variable language_">self</span>.seq.count(<span class="hljs-string">"G"</span>) + <span class="hljs-variable language_">self</span>.seq.count(<span class="hljs-string">"C"</span>)
        <span class="hljs-keyword">return</span> gc / <span class="hljs-built_in">len</span>(<span class="hljs-variable language_">self</span>.seq) * <span class="hljs-number">100</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__repr__</span>(<span class="hljs-params">self</span>):
        <span class="hljs-keyword">return</span> <span class="hljs-string">f"Sequence(<span class="hljs-subst">{self.seq_id}</span>, <span class="hljs-subst">{self.length()}</span>bp)"</span>

s = <span class="hljs-type">Sequence</span>(<span class="hljs-string">"seq1"</span>, <span class="hljs-string">"atgcgta"</span>)
<span class="hljs-built_in">print</span>(s.length())        <span class="hljs-comment"># 7</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">f"GC: <span class="hljs-subst">{s.gc_content():<span class="hljs-number">.1</span>f}</span>%"</span>)   <span class="hljs-comment"># GC: 57.1%</span>
<span class="hljs-built_in">print</span>(s)                 <span class="hljs-comment"># Sequence(seq1, 7bp)</span>
</code></pre>
<h3>3.2 Inheritance</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Protein</span>(<span class="hljs-title class_ inherited__">Sequence</span>):
    alphabet = <span class="hljs-string">"ACDEFGHIKLMNPQRSTVWY"</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, seq_id, seq</span>):
        <span class="hljs-built_in">super</span>().__init__(seq_id, seq)   <span class="hljs-comment"># Call the parent class constructor</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">molecular_weight</span>(<span class="hljs-params">self</span>):
        <span class="hljs-comment"># Simplified calculation: each amino acid is approximately 110 Da</span>
        <span class="hljs-keyword">return</span> <span class="hljs-variable language_">self</span>.length() * <span class="hljs-number">110</span>

p = Protein(<span class="hljs-string">"prot1"</span>, <span class="hljs-string">"MKWVTFISLL"</span>)
<span class="hljs-built_in">print</span>(p.molecular_weight())   <span class="hljs-comment"># 1210</span>
</code></pre>
<h3>3.3 Magic Methods</h3>
<p>Common magic methods allow objects to support built-in operations:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Vector</span>:
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, x, y</span>):
        <span class="hljs-variable language_">self</span>.x = x
        <span class="hljs-variable language_">self</span>.y = y

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__add__</span>(<span class="hljs-params">self, other</span>):      <span class="hljs-comment"># +</span>
        <span class="hljs-keyword">return</span> Vector(<span class="hljs-variable language_">self</span>.x + other.x, <span class="hljs-variable language_">self</span>.y + other.y)

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__eq__</span>(<span class="hljs-params">self, other</span>):       <span class="hljs-comment"># ==</span>
        <span class="hljs-keyword">return</span> <span class="hljs-variable language_">self</span>.x == other.x <span class="hljs-keyword">and</span> <span class="hljs-variable language_">self</span>.y == other.y

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__repr__</span>(<span class="hljs-params">self</span>):            <span class="hljs-comment"># print / repr</span>
        <span class="hljs-keyword">return</span> <span class="hljs-string">f"Vector(<span class="hljs-subst">{self.x}</span>, <span class="hljs-subst">{self.y}</span>)"</span>

v1 = Vector(<span class="hljs-number">1</span>, <span class="hljs-number">2</span>)
v2 = Vector(<span class="hljs-number">3</span>, <span class="hljs-number">4</span>)
<span class="hljs-built_in">print</span>(v1 + v2)          <span class="hljs-comment"># Vector(4, 6)</span>
<span class="hljs-built_in">print</span>(v1 == Vector(<span class="hljs-number">1</span>, <span class="hljs-number">2</span>))  <span class="hljs-comment"># True</span>
</code></pre>
<h3>3.4 Properties (property)</h3>
<p>Use @property to turn methods into attribute access, and add validation:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Person</span>:
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, name</span>):
        <span class="hljs-variable language_">self</span>._name = name

<span class="hljs-meta">    @property</span>
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">name</span>(<span class="hljs-params">self</span>):
        <span class="hljs-keyword">return</span> <span class="hljs-variable language_">self</span>._name

<span class="hljs-meta">    @name.setter</span>
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">name</span>(<span class="hljs-params">self, value</span>):
        <span class="hljs-keyword">if</span> <span class="hljs-keyword">not</span> value.strip():
            <span class="hljs-keyword">raise</span> ValueError(<span class="hljs-string">"Name cannot be empty"</span>)
        <span class="hljs-variable language_">self</span>._name = value

p = Person(<span class="hljs-string">"zorrooz"</span>)
p.name = <span class="hljs-string">"bio"</span>
<span class="hljs-built_in">print</span>(p.name)
</code></pre>
<h2>4. Exception Handling</h2>
<pre><code class="hljs language-python"><span class="hljs-keyword">try</span>:
    num = <span class="hljs-built_in">int</span>(<span class="hljs-built_in">input</span>(<span class="hljs-string">"Enter an integer: "</span>))
    result = <span class="hljs-number">100</span> / num
<span class="hljs-keyword">except</span> ValueError:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"Not an integer"</span>)
<span class="hljs-keyword">except</span> ZeroDivisionError:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"Cannot divide by zero"</span>)
<span class="hljs-keyword">else</span>:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"Result: <span class="hljs-subst">{result}</span>"</span>)
<span class="hljs-keyword">finally</span>:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"This always executes regardless of errors"</span>)
</code></pre>
<p>Custom exceptions:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">InvalidSequenceError</span>(<span class="hljs-title class_ inherited__">Exception</span>):
    <span class="hljs-keyword">pass</span>

<span class="hljs-keyword">def</span> <span class="hljs-title function_">validate</span>(<span class="hljs-params">seq</span>):
    <span class="hljs-keyword">if</span> <span class="hljs-keyword">not</span> <span class="hljs-built_in">set</span>(seq) &#x3C;= <span class="hljs-built_in">set</span>(<span class="hljs-string">"ACGT"</span>):
        <span class="hljs-keyword">raise</span> InvalidSequenceError(<span class="hljs-string">f"Contains invalid characters: <span class="hljs-subst">{seq}</span>"</span>)

<span class="hljs-keyword">try</span>:
    validate(<span class="hljs-string">"ATGXYZ"</span>)
<span class="hljs-keyword">except</span> InvalidSequenceError <span class="hljs-keyword">as</span> e:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"Validation failed: <span class="hljs-subst">{e}</span>"</span>)
</code></pre>
<h2>5. Modules and Packages</h2>
<h3>5.1 Module Imports</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> math                       <span class="hljs-comment"># Entire module</span>
<span class="hljs-keyword">from</span> math <span class="hljs-keyword">import</span> sqrt, pi         <span class="hljs-comment"># Import specific names</span>
<span class="hljs-keyword">import</span> numpy <span class="hljs-keyword">as</span> np                <span class="hljs-comment"># Alias</span>
<span class="hljs-keyword">from</span> collections <span class="hljs-keyword">import</span> Counter   <span class="hljs-comment"># Common: counting</span>

<span class="hljs-comment"># Counter example</span>
<span class="hljs-keyword">from</span> collections <span class="hljs-keyword">import</span> Counter
cnt = Counter(<span class="hljs-string">"ATGCCGA"</span>)
<span class="hljs-built_in">print</span>(cnt)          <span class="hljs-comment"># Counter({'C': 2, 'G': 2, 'A': 2, 'T': 1})</span>
<span class="hljs-built_in">print</span>(cnt.most_common(<span class="hljs-number">2</span>))
</code></pre>
<h3>5.2 Package Structure</h3>
<pre><code>myproject/
├── __init__.py        # Marks it as a package (optional in Python 3.3+)
├── utils/
│   ├── __init__.py
│   └── fasta.py       # Defines read_fasta()
└── main.py
</code></pre>
<pre><code class="hljs language-python"><span class="hljs-comment"># main.py</span>
<span class="hljs-keyword">from</span> utils.fasta <span class="hljs-keyword">import</span> read_fasta   <span class="hljs-comment"># Relative package import</span>
</code></pre>
<h3>5.3 <code>if __name__ == "__main__"</code></h3>
<p>Allows the script to be both imported and run directly:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">main</span>():
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"Running main logic"</span>)

<span class="hljs-keyword">if</span> __name__ == <span class="hljs-string">"__main__"</span>:
    main()
</code></pre>
<h2>6. Type Hints (typing)</h2>
<p>Type hints improve readability and work with IDE static checks:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">from</span> typing <span class="hljs-keyword">import</span> <span class="hljs-type">List</span>, <span class="hljs-type">Dict</span>, <span class="hljs-type">Optional</span>

<span class="hljs-keyword">def</span> <span class="hljs-title function_">count_bases</span>(<span class="hljs-params">seq: <span class="hljs-built_in">str</span></span>) -> <span class="hljs-type">Dict</span>[<span class="hljs-built_in">str</span>, <span class="hljs-built_in">int</span>]:
    <span class="hljs-string">"""Return the occurrence count of each base."""</span>
    <span class="hljs-keyword">return</span> {b: seq.count(b) <span class="hljs-keyword">for</span> b <span class="hljs-keyword">in</span> <span class="hljs-string">"ACGT"</span>}

<span class="hljs-keyword">def</span> <span class="hljs-title function_">find_motif</span>(<span class="hljs-params">seq: <span class="hljs-built_in">str</span>, motif: <span class="hljs-built_in">str</span></span>) -> <span class="hljs-type">Optional</span>[<span class="hljs-built_in">int</span>]:
    idx = seq.find(motif)
    <span class="hljs-keyword">return</span> idx <span class="hljs-keyword">if</span> idx >= <span class="hljs-number">0</span> <span class="hljs-keyword">else</span> <span class="hljs-literal">None</span>

<span class="hljs-built_in">print</span>(count_bases(<span class="hljs-string">"ATGCCGA"</span>))
</code></pre>
<h2>7. Summary</h2>
<ul>
<li><code>*args</code> / <code>**kwargs</code>, lambda, closures, and decorators are core to functional programming</li>
<li>Classes: <code>__init__</code>, inheritance, magic methods, <code>@property</code></li>
<li>Exceptions: <code>try/except/else/finally</code>, custom exceptions inherit <code>Exception</code></li>
<li>Modularity: package directories + <code>if __name__ == "__main__"</code> guard</li>
<li>Type hints make code more maintainable</li>
</ul>
<p>The next article will cover Python data processing in practice: file I/O, regular expressions, and NumPy/Pandas.</p>`,Rb=`<h1>Python 进阶：函数、类与模块</h1>
<p>掌握基础语法后，本文深入 Python 的函数式与面向对象编程，帮助你写出结构清晰、可复用的代码。</p>
<h2>1. 函数</h2>
<h3>1.1 定义与调用</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">greet</span>(<span class="hljs-params">name, greeting=<span class="hljs-string">"Hello"</span></span>):
    <span class="hljs-string">"""返回问候语。greeting 是带默认值的参数。"""</span>
    <span class="hljs-keyword">return</span> <span class="hljs-string">f"<span class="hljs-subst">{greeting}</span>, <span class="hljs-subst">{name}</span>!"</span>

<span class="hljs-built_in">print</span>(greet(<span class="hljs-string">"zorrooz"</span>))            <span class="hljs-comment"># Hello, zorrooz!</span>
<span class="hljs-built_in">print</span>(greet(<span class="hljs-string">"zorrooz"</span>, <span class="hljs-string">"Hi"</span>))      <span class="hljs-comment"># Hi, zorrooz!</span>
</code></pre>
<h3>1.2 参数传递</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">show</span>(<span class="hljs-params">a, b, *args, **kwargs</span>):
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"a=<span class="hljs-subst">{a}</span>, b=<span class="hljs-subst">{b}</span>"</span>)
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"位置参数: <span class="hljs-subst">{args}</span>"</span>)      <span class="hljs-comment"># 元组</span>
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"关键字参数: <span class="hljs-subst">{kwargs}</span>"</span>)  <span class="hljs-comment"># 字典</span>

show(<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>, x=<span class="hljs-number">5</span>, y=<span class="hljs-number">6</span>)
<span class="hljs-comment"># a=1, b=2</span>
<span class="hljs-comment"># 位置参数: (3, 4)</span>
<span class="hljs-comment"># 关键字参数: {'x': 5, 'y': 6}</span>
</code></pre>
<ul>
<li><code>*args</code>：收集多余的位置参数为元组</li>
<li><code>**kwargs</code>：收集多余的关键字参数为字典</li>
</ul>
<h3>1.3 解包传参</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">add</span>(<span class="hljs-params">a, b, c</span>):
    <span class="hljs-keyword">return</span> a + b + c

nums = [<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>]
<span class="hljs-built_in">print</span>(add(*nums))            <span class="hljs-comment"># 6，列表解包</span>

data = {<span class="hljs-string">"a"</span>: <span class="hljs-number">1</span>, <span class="hljs-string">"b"</span>: <span class="hljs-number">2</span>, <span class="hljs-string">"c"</span>: <span class="hljs-number">3</span>}
<span class="hljs-built_in">print</span>(add(**data))           <span class="hljs-comment"># 6，字典解包</span>
</code></pre>
<h3>1.4 lambda 匿名函数</h3>
<pre><code class="hljs language-python">square = <span class="hljs-keyword">lambda</span> x: x ** <span class="hljs-number">2</span>
<span class="hljs-built_in">print</span>(square(<span class="hljs-number">5</span>))             <span class="hljs-comment"># 25</span>

<span class="hljs-comment"># 与 sorted / map / filter 配合</span>
words = [<span class="hljs-string">"banana"</span>, <span class="hljs-string">"apple"</span>, <span class="hljs-string">"cherry"</span>]
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">sorted</span>(words, key=<span class="hljs-keyword">lambda</span> w: <span class="hljs-built_in">len</span>(w)))
<span class="hljs-comment"># ['apple', 'banana', 'cherry']</span>

nums = [<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>, <span class="hljs-number">5</span>]
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">list</span>(<span class="hljs-built_in">map</span>(<span class="hljs-keyword">lambda</span> x: x * <span class="hljs-number">2</span>, nums)))       <span class="hljs-comment"># [2, 4, 6, 8, 10]</span>
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">list</span>(<span class="hljs-built_in">filter</span>(<span class="hljs-keyword">lambda</span> x: x % <span class="hljs-number">2</span> == <span class="hljs-number">0</span>, nums)))  <span class="hljs-comment"># [2, 4]</span>
</code></pre>
<h3>1.5 闭包</h3>
<p>函数内部定义函数并引用外部变量：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">make_counter</span>():
    count = <span class="hljs-number">0</span>
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">counter</span>():
        <span class="hljs-keyword">nonlocal</span> count
        count += <span class="hljs-number">1</span>
        <span class="hljs-keyword">return</span> count
    <span class="hljs-keyword">return</span> counter

c = make_counter()
<span class="hljs-built_in">print</span>(c())   <span class="hljs-comment"># 1</span>
<span class="hljs-built_in">print</span>(c())   <span class="hljs-comment"># 2</span>
</code></pre>
<p><code>nonlocal</code> 用于修改外层函数中的变量。</p>
<h2>2. 装饰器</h2>
<p>装饰器是"接收函数、返回函数"的高阶函数，用于在不修改原函数的情况下增强其行为：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> time

<span class="hljs-keyword">def</span> <span class="hljs-title function_">timer</span>(<span class="hljs-params">func</span>):
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">wrapper</span>(<span class="hljs-params">*args, **kwargs</span>):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{func.__name__}</span> 耗时 <span class="hljs-subst">{elapsed:<span class="hljs-number">.4</span>f}</span>s"</span>)
        <span class="hljs-keyword">return</span> result
    <span class="hljs-keyword">return</span> wrapper

<span class="hljs-meta">@timer</span>
<span class="hljs-keyword">def</span> <span class="hljs-title function_">slow_sum</span>(<span class="hljs-params">n</span>):
    <span class="hljs-keyword">return</span> <span class="hljs-built_in">sum</span>(<span class="hljs-built_in">range</span>(n))

<span class="hljs-built_in">print</span>(slow_sum(<span class="hljs-number">10_000_000</span>))
<span class="hljs-comment"># slow_sum 耗时 0.35xx s</span>
</code></pre>
<p>带参数的装饰器：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">repeat</span>(<span class="hljs-params">times</span>):
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">decorator</span>(<span class="hljs-params">func</span>):
        <span class="hljs-keyword">def</span> <span class="hljs-title function_">wrapper</span>(<span class="hljs-params">*args, **kwargs</span>):
            <span class="hljs-keyword">for</span> _ <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(times):
                func(*args, **kwargs)
        <span class="hljs-keyword">return</span> wrapper
    <span class="hljs-keyword">return</span> decorator

<span class="hljs-meta">@repeat(<span class="hljs-params"><span class="hljs-number">3</span></span>)</span>
<span class="hljs-keyword">def</span> <span class="hljs-title function_">hello</span>():
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"Hi!"</span>)

hello()
<span class="hljs-comment"># Hi! Hi! Hi!</span>
</code></pre>
<h2>3. 类与面向对象</h2>
<h3>3.1 基本定义</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Sequence</span>:
    <span class="hljs-string">"""表示一条生物序列。"""</span>

    <span class="hljs-comment"># 类属性：所有实例共享</span>
    alphabet = <span class="hljs-string">"ACGT"</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, seq_id, seq</span>):
        <span class="hljs-string">"""构造方法：初始化实例属性。"""</span>
        <span class="hljs-variable language_">self</span>.seq_id = seq_id
        <span class="hljs-variable language_">self</span>.seq = seq.upper()

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">length</span>(<span class="hljs-params">self</span>):
        <span class="hljs-keyword">return</span> <span class="hljs-built_in">len</span>(<span class="hljs-variable language_">self</span>.seq)

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">gc_content</span>(<span class="hljs-params">self</span>):
        gc = <span class="hljs-variable language_">self</span>.seq.count(<span class="hljs-string">"G"</span>) + <span class="hljs-variable language_">self</span>.seq.count(<span class="hljs-string">"C"</span>)
        <span class="hljs-keyword">return</span> gc / <span class="hljs-built_in">len</span>(<span class="hljs-variable language_">self</span>.seq) * <span class="hljs-number">100</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__repr__</span>(<span class="hljs-params">self</span>):
        <span class="hljs-keyword">return</span> <span class="hljs-string">f"Sequence(<span class="hljs-subst">{self.seq_id}</span>, <span class="hljs-subst">{self.length()}</span>bp)"</span>

s = <span class="hljs-type">Sequence</span>(<span class="hljs-string">"seq1"</span>, <span class="hljs-string">"atgcgta"</span>)
<span class="hljs-built_in">print</span>(s.length())        <span class="hljs-comment"># 7</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">f"GC: <span class="hljs-subst">{s.gc_content():<span class="hljs-number">.1</span>f}</span>%"</span>)   <span class="hljs-comment"># GC: 57.1%</span>
<span class="hljs-built_in">print</span>(s)                 <span class="hljs-comment"># Sequence(seq1, 7bp)</span>
</code></pre>
<h3>3.2 继承</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Protein</span>(<span class="hljs-title class_ inherited__">Sequence</span>):
    alphabet = <span class="hljs-string">"ACDEFGHIKLMNPQRSTVWY"</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, seq_id, seq</span>):
        <span class="hljs-built_in">super</span>().__init__(seq_id, seq)   <span class="hljs-comment"># 调用父类构造方法</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">molecular_weight</span>(<span class="hljs-params">self</span>):
        <span class="hljs-comment"># 简化计算：每个氨基酸约 110 Da</span>
        <span class="hljs-keyword">return</span> <span class="hljs-variable language_">self</span>.length() * <span class="hljs-number">110</span>

p = Protein(<span class="hljs-string">"prot1"</span>, <span class="hljs-string">"MKWVTFISLL"</span>)
<span class="hljs-built_in">print</span>(p.molecular_weight())   <span class="hljs-comment"># 1210</span>
</code></pre>
<h3>3.3 魔术方法</h3>
<p>常用魔术方法让对象支持内置操作：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Vector</span>:
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, x, y</span>):
        <span class="hljs-variable language_">self</span>.x = x
        <span class="hljs-variable language_">self</span>.y = y

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__add__</span>(<span class="hljs-params">self, other</span>):      <span class="hljs-comment"># +</span>
        <span class="hljs-keyword">return</span> Vector(<span class="hljs-variable language_">self</span>.x + other.x, <span class="hljs-variable language_">self</span>.y + other.y)

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__eq__</span>(<span class="hljs-params">self, other</span>):       <span class="hljs-comment"># ==</span>
        <span class="hljs-keyword">return</span> <span class="hljs-variable language_">self</span>.x == other.x <span class="hljs-keyword">and</span> <span class="hljs-variable language_">self</span>.y == other.y

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__repr__</span>(<span class="hljs-params">self</span>):            <span class="hljs-comment"># print / repr</span>
        <span class="hljs-keyword">return</span> <span class="hljs-string">f"Vector(<span class="hljs-subst">{self.x}</span>, <span class="hljs-subst">{self.y}</span>)"</span>

v1 = Vector(<span class="hljs-number">1</span>, <span class="hljs-number">2</span>)
v2 = Vector(<span class="hljs-number">3</span>, <span class="hljs-number">4</span>)
<span class="hljs-built_in">print</span>(v1 + v2)          <span class="hljs-comment"># Vector(4, 6)</span>
<span class="hljs-built_in">print</span>(v1 == Vector(<span class="hljs-number">1</span>, <span class="hljs-number">2</span>))  <span class="hljs-comment"># True</span>
</code></pre>
<h3>3.4 属性（property）</h3>
<p>用 <code>@property</code> 把方法变成属性访问，可加入校验：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Person</span>:
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, name</span>):
        <span class="hljs-variable language_">self</span>._name = name

<span class="hljs-meta">    @property</span>
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">name</span>(<span class="hljs-params">self</span>):
        <span class="hljs-keyword">return</span> <span class="hljs-variable language_">self</span>._name

<span class="hljs-meta">    @name.setter</span>
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">name</span>(<span class="hljs-params">self, value</span>):
        <span class="hljs-keyword">if</span> <span class="hljs-keyword">not</span> value.strip():
            <span class="hljs-keyword">raise</span> ValueError(<span class="hljs-string">"名字不能为空"</span>)
        <span class="hljs-variable language_">self</span>._name = value

p = Person(<span class="hljs-string">"zorrooz"</span>)
p.name = <span class="hljs-string">"bio"</span>
<span class="hljs-built_in">print</span>(p.name)
</code></pre>
<h2>4. 异常处理</h2>
<pre><code class="hljs language-python"><span class="hljs-keyword">try</span>:
    num = <span class="hljs-built_in">int</span>(<span class="hljs-built_in">input</span>(<span class="hljs-string">"输入一个整数："</span>))
    result = <span class="hljs-number">100</span> / num
<span class="hljs-keyword">except</span> ValueError:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"输入的不是整数"</span>)
<span class="hljs-keyword">except</span> ZeroDivisionError:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"不能除以零"</span>)
<span class="hljs-keyword">else</span>:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"结果: <span class="hljs-subst">{result}</span>"</span>)
<span class="hljs-keyword">finally</span>:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"无论是否出错都会执行"</span>)
</code></pre>
<p>自定义异常：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">InvalidSequenceError</span>(<span class="hljs-title class_ inherited__">Exception</span>):
    <span class="hljs-keyword">pass</span>

<span class="hljs-keyword">def</span> <span class="hljs-title function_">validate</span>(<span class="hljs-params">seq</span>):
    <span class="hljs-keyword">if</span> <span class="hljs-keyword">not</span> <span class="hljs-built_in">set</span>(seq) &#x3C;= <span class="hljs-built_in">set</span>(<span class="hljs-string">"ACGT"</span>):
        <span class="hljs-keyword">raise</span> InvalidSequenceError(<span class="hljs-string">f"包含非法字符: <span class="hljs-subst">{seq}</span>"</span>)

<span class="hljs-keyword">try</span>:
    validate(<span class="hljs-string">"ATGXYZ"</span>)
<span class="hljs-keyword">except</span> InvalidSequenceError <span class="hljs-keyword">as</span> e:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"校验失败: <span class="hljs-subst">{e}</span>"</span>)
</code></pre>
<h2>5. 模块与包</h2>
<h3>5.1 模块导入</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> math                       <span class="hljs-comment"># 整个模块</span>
<span class="hljs-keyword">from</span> math <span class="hljs-keyword">import</span> sqrt, pi         <span class="hljs-comment"># 导入指定名字</span>
<span class="hljs-keyword">import</span> numpy <span class="hljs-keyword">as</span> np                <span class="hljs-comment"># 别名</span>
<span class="hljs-keyword">from</span> collections <span class="hljs-keyword">import</span> Counter   <span class="hljs-comment"># 常用：计数</span>

<span class="hljs-comment"># Counter 示例</span>
<span class="hljs-keyword">from</span> collections <span class="hljs-keyword">import</span> Counter
cnt = Counter(<span class="hljs-string">"ATGCCGA"</span>)
<span class="hljs-built_in">print</span>(cnt)          <span class="hljs-comment"># Counter({'C': 2, 'G': 2, 'A': 2, 'T': 1})</span>
<span class="hljs-built_in">print</span>(cnt.most_common(<span class="hljs-number">2</span>))
</code></pre>
<h3>5.2 包结构</h3>
<pre><code>myproject/
├── __init__.py        # 标记为包（Python 3.3+ 可省略）
├── utils/
│   ├── __init__.py
│   └── fasta.py       # 定义 read_fasta()
└── main.py
</code></pre>
<pre><code class="hljs language-python"><span class="hljs-comment"># main.py</span>
<span class="hljs-keyword">from</span> utils.fasta <span class="hljs-keyword">import</span> read_fasta   <span class="hljs-comment"># 相对包路径导入</span>
</code></pre>
<h3>5.3 <code>if __name__ == "__main__"</code></h3>
<p>让脚本既可被导入也可直接运行：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">main</span>():
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"运行主逻辑"</span>)

<span class="hljs-keyword">if</span> __name__ == <span class="hljs-string">"__main__"</span>:
    main()
</code></pre>
<h2>6. 类型提示（typing）</h2>
<p>类型提示提高可读性，配合 IDE 静态检查：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">from</span> typing <span class="hljs-keyword">import</span> <span class="hljs-type">List</span>, <span class="hljs-type">Dict</span>, <span class="hljs-type">Optional</span>

<span class="hljs-keyword">def</span> <span class="hljs-title function_">count_bases</span>(<span class="hljs-params">seq: <span class="hljs-built_in">str</span></span>) -> <span class="hljs-type">Dict</span>[<span class="hljs-built_in">str</span>, <span class="hljs-built_in">int</span>]:
    <span class="hljs-string">"""返回每种碱基的出现次数。"""</span>
    <span class="hljs-keyword">return</span> {b: seq.count(b) <span class="hljs-keyword">for</span> b <span class="hljs-keyword">in</span> <span class="hljs-string">"ACGT"</span>}

<span class="hljs-keyword">def</span> <span class="hljs-title function_">find_motif</span>(<span class="hljs-params">seq: <span class="hljs-built_in">str</span>, motif: <span class="hljs-built_in">str</span></span>) -> <span class="hljs-type">Optional</span>[<span class="hljs-built_in">int</span>]:
    idx = seq.find(motif)
    <span class="hljs-keyword">return</span> idx <span class="hljs-keyword">if</span> idx >= <span class="hljs-number">0</span> <span class="hljs-keyword">else</span> <span class="hljs-literal">None</span>

<span class="hljs-built_in">print</span>(count_bases(<span class="hljs-string">"ATGCCGA"</span>))
</code></pre>
<h2>7. 小结</h2>
<ul>
<li><code>*args</code> / <code>**kwargs</code>、lambda、闭包、装饰器是函数式编程核心</li>
<li>类：<code>__init__</code>、继承、魔术方法、<code>@property</code></li>
<li>异常：<code>try/except/else/finally</code>，自定义异常继承 <code>Exception</code></li>
<li>模块化：包目录 + <code>if __name__ == "__main__"</code> 守卫</li>
<li>类型提示让代码更可维护</li>
</ul>
<p>下一篇将介绍 Python 数据处理实战：文件 IO、正则表达式与 NumPy/Pandas。</p>`,Lb=`<h1>Introduction to Python Programming: Environment, Syntax, and Data Types</h1>
<p>Python is an interpreted language with clean syntax and a rich ecosystem, widely used in scientific computing, data processing, and script automation. This article starts from scratch, covering environment setup, basic syntax, built-in data types, and control flow. All examples are directly runnable.</p>
<h2>1. Environment Setup</h2>
<h3>1.1 Installing Python</h3>
<p>Visit <a href="https://www.python.org/">python.org</a> to download the installer for your platform. When installing, be sure to check <strong>Add Python to PATH</strong>.</p>
<p>Verify the installation:</p>
<pre><code class="hljs language-bash">python --version
<span class="hljs-comment"># Python 3.12.x</span>
</code></pre>
<blockquote>
<p>On Linux / WSL, the <code>python3</code> command is typically used; on Windows, it is <code>python</code>.</p>
</blockquote>
<h3>1.2 Virtual Environment (Highly Recommended)</h3>
<p>A virtual environment isolates dependencies for each project, avoiding pollution of the global environment:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Create a virtual environment</span>
python -m venv .venv

<span class="hljs-comment"># Activate (Windows)</span>
.venv\\Scripts\\activate
<span class="hljs-comment"># Activate (Linux / macOS / WSL)</span>
<span class="hljs-built_in">source</span> .venv/bin/activate
</code></pre>
<p>After activation, <code>(.venv)</code> appears before the prompt. Install dependencies:</p>
<pre><code class="hljs language-bash">pip install numpy pandas
</code></pre>
<p>Export the dependency list to a file for easy reproduction:</p>
<pre><code class="hljs language-bash">pip freeze > requirements.txt
</code></pre>
<h3>1.3 Interactive Interpreter and Scripts</h3>
<pre><code class="hljs language-bash">python            <span class="hljs-comment"># Enter REPL interactive mode</span>
python script.py  <span class="hljs-comment"># Run a script file</span>
</code></pre>
<p>In the REPL, you can directly enter expressions and get immediate results, which is great for quickly validating ideas.</p>
<h2>2. First Program</h2>
<pre><code class="hljs language-python"><span class="hljs-built_in">print</span>(<span class="hljs-string">"Hello, Python!"</span>)
</code></pre>
<p>Save it as <code>hello.py</code> and run:</p>
<pre><code class="hljs language-bash">python hello.py
<span class="hljs-comment"># Hello, Python!</span>
</code></pre>
<h2>3. Basic Syntax</h2>
<h3>3.1 Comments</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Single-line comment</span>

<span class="hljs-string">"""
Multi-line comment (actually a multi-line string,
often used as documentation)
"""</span>
</code></pre>
<h3>3.2 Variables and Assignment</h3>
<p>Python is a dynamically typed language; variables do not require type declarations:</p>
<pre><code class="hljs language-python">name = <span class="hljs-string">"zorrooz"</span>     <span class="hljs-comment"># string</span>
age = <span class="hljs-number">25</span>             <span class="hljs-comment"># integer</span>
height = <span class="hljs-number">1.78</span>        <span class="hljs-comment"># float</span>
is_student = <span class="hljs-literal">True</span>    <span class="hljs-comment"># boolean</span>

<span class="hljs-comment"># Assign multiple variables at once</span>
a, b, c = <span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>

<span class="hljs-comment"># Swap two variables</span>
a, b = b, a
</code></pre>
<p>Variable naming follows PEP 8: lowercase letters + underscores (<code>snake_case</code>).</p>
<h3>3.3 Input and Output</h3>
<pre><code class="hljs language-python">name = <span class="hljs-built_in">input</span>(<span class="hljs-string">"请输入你的名字："</span>)
<span class="hljs-built_in">print</span>(<span class="hljs-string">"你好，"</span>, name)

<span class="hljs-comment"># f-string formatting (Python 3.6+, recommended)</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">f"你好，<span class="hljs-subst">{name}</span>，今年 <span class="hljs-subst">{age}</span> 岁"</span>)

<span class="hljs-comment"># format method</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">"你好，{}，今年 {} 岁"</span>.<span class="hljs-built_in">format</span>(name, age))
</code></pre>
<h2>4. Built-in Data Types</h2>
<h3>4.1 Numbers</h3>
<pre><code class="hljs language-python">x = <span class="hljs-number">10</span>          <span class="hljs-comment"># int</span>
y = <span class="hljs-number">3.14</span>        <span class="hljs-comment"># float</span>
z = <span class="hljs-number">2</span> + <span class="hljs-number">3j</span>      <span class="hljs-comment"># complex</span>

<span class="hljs-comment"># Common operations</span>
<span class="hljs-built_in">print</span>(<span class="hljs-number">7</span> // <span class="hljs-number">2</span>)   <span class="hljs-comment"># integer division, result 3</span>
<span class="hljs-built_in">print</span>(<span class="hljs-number">7</span> % <span class="hljs-number">2</span>)    <span class="hljs-comment"># modulo, result 1</span>
<span class="hljs-built_in">print</span>(<span class="hljs-number">2</span> ** <span class="hljs-number">10</span>)  <span class="hljs-comment"># exponentiation, result 1024</span>
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">round</span>(<span class="hljs-number">3.14159</span>, <span class="hljs-number">2</span>))  <span class="hljs-comment"># round to 2 decimal places</span>
</code></pre>
<h3>4.2 Strings</h3>
<pre><code class="hljs language-python">s = <span class="hljs-string">"hello"</span>
t = <span class="hljs-string">'world'</span>

<span class="hljs-comment"># Concatenation and repetition</span>
<span class="hljs-built_in">print</span>(s + <span class="hljs-string">" "</span> + t)   <span class="hljs-comment"># hello world</span>
<span class="hljs-built_in">print</span>(s * <span class="hljs-number">3</span>)         <span class="hljs-comment"># hellohellohello</span>

<span class="hljs-comment"># Indexing and slicing (left-inclusive, right-exclusive)</span>
<span class="hljs-built_in">print</span>(s[<span class="hljs-number">0</span>])          <span class="hljs-comment"># h</span>
<span class="hljs-built_in">print</span>(s[-<span class="hljs-number">1</span>])         <span class="hljs-comment"># o</span>
<span class="hljs-built_in">print</span>(s[<span class="hljs-number">1</span>:<span class="hljs-number">3</span>])        <span class="hljs-comment"># el</span>

<span class="hljs-comment"># Common methods</span>
<span class="hljs-built_in">print</span>(s.upper())          <span class="hljs-comment"># HELLO</span>
<span class="hljs-built_in">print</span>(s.replace(<span class="hljs-string">"l"</span>, <span class="hljs-string">"L"</span>))  <span class="hljs-comment"># heLLo</span>
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">len</span>(s))             <span class="hljs-comment"># 5</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">", "</span>.join([<span class="hljs-string">"a"</span>, <span class="hljs-string">"b"</span>, <span class="hljs-string">"c"</span>]))  <span class="hljs-comment"># a, b, c</span>
</code></pre>
<h3>4.3 Lists (list)</h3>
<p>Lists are ordered, mutable containers:</p>
<pre><code class="hljs language-python">fruits = [<span class="hljs-string">"apple"</span>, <span class="hljs-string">"banana"</span>, <span class="hljs-string">"cherry"</span>]

fruits.append(<span class="hljs-string">"orange"</span>)      <span class="hljs-comment"># append at the end</span>
fruits.insert(<span class="hljs-number">0</span>, <span class="hljs-string">"grape"</span>)    <span class="hljs-comment"># insert at a specific position</span>
fruits.remove(<span class="hljs-string">"apple"</span>)       <span class="hljs-comment"># remove by value</span>
popped = fruits.pop()        <span class="hljs-comment"># pop the last element</span>

<span class="hljs-built_in">print</span>(fruits[<span class="hljs-number">0</span>])             <span class="hljs-comment"># indexing</span>
<span class="hljs-built_in">print</span>(fruits[-<span class="hljs-number">1</span>])            <span class="hljs-comment"># last element</span>
<span class="hljs-built_in">print</span>(fruits[<span class="hljs-number">1</span>:<span class="hljs-number">3</span>])           <span class="hljs-comment"># slicing</span>

<span class="hljs-comment"># List comprehension (very commonly used)</span>
squares = [x ** <span class="hljs-number">2</span> <span class="hljs-keyword">for</span> x <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">10</span>)]
<span class="hljs-built_in">print</span>(squares)  <span class="hljs-comment"># [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]</span>
</code></pre>
<h3>4.4 Tuples (tuple)</h3>
<p>Tuples are immutable lists, often used for returning multiple values:</p>
<pre><code class="hljs language-python">point = (<span class="hljs-number">3</span>, <span class="hljs-number">4</span>)
x, y = point          <span class="hljs-comment"># unpacking</span>
<span class="hljs-built_in">print</span>(x, y)           <span class="hljs-comment"># 3 4</span>

<span class="hljs-keyword">def</span> <span class="hljs-title function_">minmax</span>(<span class="hljs-params">nums</span>):
    <span class="hljs-keyword">return</span> <span class="hljs-built_in">min</span>(nums), <span class="hljs-built_in">max</span>(nums)

lo, hi = minmax([<span class="hljs-number">3</span>, <span class="hljs-number">1</span>, <span class="hljs-number">4</span>, <span class="hljs-number">1</span>, <span class="hljs-number">5</span>])
<span class="hljs-built_in">print</span>(lo, hi)         <span class="hljs-comment"># 1 5</span>
</code></pre>
<h3>4.5 Dictionaries (dict)</h3>
<p>Dictionaries are key-value containers with efficient lookup:</p>
<pre><code class="hljs language-python">person = {<span class="hljs-string">"name"</span>: <span class="hljs-string">"zorrooz"</span>, <span class="hljs-string">"age"</span>: <span class="hljs-number">25</span>}

<span class="hljs-built_in">print</span>(person[<span class="hljs-string">"name"</span>])            <span class="hljs-comment"># access value</span>
<span class="hljs-built_in">print</span>(person.get(<span class="hljs-string">"city"</span>, <span class="hljs-string">"未知"</span>))  <span class="hljs-comment"># safe access with default value</span>
person[<span class="hljs-string">"city"</span>] = <span class="hljs-string">"Lanzhou"</span>       <span class="hljs-comment"># add/update</span>

<span class="hljs-keyword">for</span> key, value <span class="hljs-keyword">in</span> person.items():
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{key}</span>: <span class="hljs-subst">{value}</span>"</span>)

<span class="hljs-comment"># Dictionary comprehension</span>
squares = {x: x ** <span class="hljs-number">2</span> <span class="hljs-keyword">for</span> x <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">5</span>)}
<span class="hljs-built_in">print</span>(squares)  <span class="hljs-comment"># {0: 0, 1: 1, 2: 4, 3: 9, 4: 16}</span>
</code></pre>
<h3>4.6 Sets (set)</h3>
<p>Sets are unordered, deduplicating containers:</p>
<pre><code class="hljs language-python">a = {<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">3</span>, <span class="hljs-number">3</span>}
<span class="hljs-built_in">print</span>(a)              <span class="hljs-comment"># {1, 2, 3}, automatically deduplicated</span>

b = {<span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>}
<span class="hljs-built_in">print</span>(a &#x26; b)          <span class="hljs-comment"># intersection {2, 3}</span>
<span class="hljs-built_in">print</span>(a | b)          <span class="hljs-comment"># union {1, 2, 3, 4}</span>
<span class="hljs-built_in">print</span>(a - b)          <span class="hljs-comment"># difference {1}</span>
</code></pre>
<h2>5. Control Flow</h2>
<h3>5.1 Conditional Statements</h3>
<pre><code class="hljs language-python">score = <span class="hljs-number">85</span>

<span class="hljs-keyword">if</span> score >= <span class="hljs-number">90</span>:
    grade = <span class="hljs-string">"A"</span>
<span class="hljs-keyword">elif</span> score >= <span class="hljs-number">80</span>:
    grade = <span class="hljs-string">"B"</span>
<span class="hljs-keyword">elif</span> score >= <span class="hljs-number">70</span>:
    grade = <span class="hljs-string">"C"</span>
<span class="hljs-keyword">else</span>:
    grade = <span class="hljs-string">"D"</span>

<span class="hljs-built_in">print</span>(<span class="hljs-string">f"等级：<span class="hljs-subst">{grade}</span>"</span>)

<span class="hljs-comment"># Ternary expression</span>
status = <span class="hljs-string">"pass"</span> <span class="hljs-keyword">if</span> score >= <span class="hljs-number">60</span> <span class="hljs-keyword">else</span> <span class="hljs-string">"fail"</span>
</code></pre>
<h3>5.2 for Loops</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Iterate over iterable objects</span>
<span class="hljs-keyword">for</span> fruit <span class="hljs-keyword">in</span> [<span class="hljs-string">"apple"</span>, <span class="hljs-string">"banana"</span>]:
    <span class="hljs-built_in">print</span>(fruit)

<span class="hljs-comment"># Generate number sequences with range</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">5</span>):        <span class="hljs-comment"># 0 to 4</span>
    <span class="hljs-built_in">print</span>(i)

<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">1</span>, <span class="hljs-number">10</span>, <span class="hljs-number">2</span>):  <span class="hljs-comment"># 1, 3, 5, 7, 9</span>
    <span class="hljs-built_in">print</span>(i)

<span class="hljs-comment"># Use enumerate to get index and value simultaneously</span>
<span class="hljs-keyword">for</span> idx, fruit <span class="hljs-keyword">in</span> <span class="hljs-built_in">enumerate</span>([<span class="hljs-string">"a"</span>, <span class="hljs-string">"b"</span>]):
    <span class="hljs-built_in">print</span>(idx, fruit)
</code></pre>
<h3>5.3 while Loops</h3>
<pre><code class="hljs language-python">n = <span class="hljs-number">0</span>
<span class="hljs-keyword">while</span> n &#x3C; <span class="hljs-number">5</span>:
    <span class="hljs-built_in">print</span>(n)
    n += <span class="hljs-number">1</span>          <span class="hljs-comment"># Note: forgetting to increment causes an infinite loop</span>
</code></pre>
<h3>5.4 break / continue / else</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># break: terminate the loop early</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">10</span>):
    <span class="hljs-keyword">if</span> i == <span class="hljs-number">3</span>:
        <span class="hljs-keyword">break</span>
    <span class="hljs-built_in">print</span>(i)        <span class="hljs-comment"># outputs 0 1 2</span>

<span class="hljs-comment"># continue: skip the current iteration</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">5</span>):
    <span class="hljs-keyword">if</span> i == <span class="hljs-number">2</span>:
        <span class="hljs-keyword">continue</span>
    <span class="hljs-built_in">print</span>(i)        <span class="hljs-comment"># outputs 0 1 3 4</span>

<span class="hljs-comment"># for-else: executes when the loop completes normally (not broken)</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">3</span>):
    <span class="hljs-keyword">if</span> i == <span class="hljs-number">99</span>:
        <span class="hljs-keyword">break</span>
<span class="hljs-keyword">else</span>:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"循环正常结束"</span>)
</code></pre>
<h2>6. Comprehensive Exercise</h2>
<p>Count the occurrences of each word in a text:</p>
<pre><code class="hljs language-python">text = <span class="hljs-string">"python is fun and python is powerful"</span>

words = text.split()
counter = {}

<span class="hljs-keyword">for</span> word <span class="hljs-keyword">in</span> words:
    counter[word] = counter.get(word, <span class="hljs-number">0</span>) + <span class="hljs-number">1</span>

<span class="hljs-keyword">for</span> word, count <span class="hljs-keyword">in</span> <span class="hljs-built_in">sorted</span>(counter.items(), key=<span class="hljs-keyword">lambda</span> x: -x[<span class="hljs-number">1</span>]):
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{word}</span>: <span class="hljs-subst">{count}</span>"</span>)
</code></pre>
<p>Output:</p>
<pre><code>python: 2
is: 2
fun: 1
and: 1
powerful: 1
</code></pre>
<h2>7. Summary</h2>
<ul>
<li>Use <code>venv</code> to isolate project dependencies, use <code>pip freeze</code> to record dependencies</li>
<li>Built-in types: <code>int</code> / <code>float</code> / <code>str</code> / <code>list</code> / <code>tuple</code> / <code>dict</code> / <code>set</code></li>
<li>List and dictionary comprehensions and f-strings are high-frequency syntax in daily use</li>
<li>Control flow: <code>if/elif/else</code>, <code>for</code>, <code>while</code>, <code>break</code> / <code>continue</code></li>
</ul>
<p>The next article will introduce functions, classes, and modules, moving into real engineering-style programming.</p>`,Mb=`<h1>Python 编程入门：环境、语法与数据类型</h1>
<p>Python 是一门语法简洁、生态丰富的解释型语言，广泛应用于科学计算、数据处理与脚本自动化。本文从零开始，覆盖环境搭建、基础语法、内置数据类型与流程控制，所有示例均可直接运行。</p>
<h2>1. 环境搭建</h2>
<h3>1.1 安装 Python</h3>
<p>访问 <a href="https://www.python.org/">python.org</a> 下载对应平台的安装包。安装时务必勾选 <strong>Add Python to PATH</strong>。</p>
<p>验证安装：</p>
<pre><code class="hljs language-bash">python --version
<span class="hljs-comment"># Python 3.12.x</span>
</code></pre>
<blockquote>
<p>在 Linux / WSL 下通常使用 <code>python3</code> 命令，Windows 下为 <code>python</code>。</p>
</blockquote>
<h3>1.2 虚拟环境（强烈推荐）</h3>
<p>虚拟环境为每个项目隔离依赖，避免污染全局环境：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 创建虚拟环境</span>
python -m venv .venv

<span class="hljs-comment"># 激活（Windows）</span>
.venv\\Scripts\\activate
<span class="hljs-comment"># 激活（Linux / macOS / WSL）</span>
<span class="hljs-built_in">source</span> .venv/bin/activate
</code></pre>
<p>激活后提示符前会出现 <code>(.venv)</code>。安装依赖：</p>
<pre><code class="hljs language-bash">pip install numpy pandas
</code></pre>
<p>将依赖清单导出到文件，方便复现：</p>
<pre><code class="hljs language-bash">pip freeze > requirements.txt
</code></pre>
<h3>1.3 交互式解释器与脚本</h3>
<pre><code class="hljs language-bash">python            <span class="hljs-comment"># 进入 REPL 交互模式</span>
python script.py  <span class="hljs-comment"># 运行脚本文件</span>
</code></pre>
<p>在 REPL 中可直接输入表达式并立即得到结果，适合快速验证想法。</p>
<h2>2. 第一个程序</h2>
<pre><code class="hljs language-python"><span class="hljs-built_in">print</span>(<span class="hljs-string">"Hello, Python!"</span>)
</code></pre>
<p>保存为 <code>hello.py</code> 后运行：</p>
<pre><code class="hljs language-bash">python hello.py
<span class="hljs-comment"># Hello, Python!</span>
</code></pre>
<h2>3. 基础语法</h2>
<h3>3.1 注释</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 单行注释</span>

<span class="hljs-string">"""
多行注释（实际上是多行字符串，
但常被用作文档说明）
"""</span>
</code></pre>
<h3>3.2 变量与赋值</h3>
<p>Python 是动态类型语言，变量无需声明类型：</p>
<pre><code class="hljs language-python">name = <span class="hljs-string">"zorrooz"</span>     <span class="hljs-comment"># 字符串</span>
age = <span class="hljs-number">25</span>             <span class="hljs-comment"># 整数</span>
height = <span class="hljs-number">1.78</span>        <span class="hljs-comment"># 浮点数</span>
is_student = <span class="hljs-literal">True</span>    <span class="hljs-comment"># 布尔值</span>

<span class="hljs-comment"># 同时赋值多个变量</span>
a, b, c = <span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>

<span class="hljs-comment"># 交换两个变量的值</span>
a, b = b, a
</code></pre>
<p>变量命名遵循 PEP 8：小写字母 + 下划线（<code>snake_case</code>）。</p>
<h3>3.3 输入输出</h3>
<pre><code class="hljs language-python">name = <span class="hljs-built_in">input</span>(<span class="hljs-string">"请输入你的名字："</span>)
<span class="hljs-built_in">print</span>(<span class="hljs-string">"你好，"</span>, name)

<span class="hljs-comment"># f-string 格式化（Python 3.6+，最推荐）</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">f"你好，<span class="hljs-subst">{name}</span>，今年 <span class="hljs-subst">{age}</span> 岁"</span>)

<span class="hljs-comment"># format 方法</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">"你好，{}，今年 {} 岁"</span>.<span class="hljs-built_in">format</span>(name, age))
</code></pre>
<h2>4. 内置数据类型</h2>
<h3>4.1 数字</h3>
<pre><code class="hljs language-python">x = <span class="hljs-number">10</span>          <span class="hljs-comment"># int</span>
y = <span class="hljs-number">3.14</span>        <span class="hljs-comment"># float</span>
z = <span class="hljs-number">2</span> + <span class="hljs-number">3j</span>      <span class="hljs-comment"># complex</span>

<span class="hljs-comment"># 常用运算</span>
<span class="hljs-built_in">print</span>(<span class="hljs-number">7</span> // <span class="hljs-number">2</span>)   <span class="hljs-comment"># 整除，得 3</span>
<span class="hljs-built_in">print</span>(<span class="hljs-number">7</span> % <span class="hljs-number">2</span>)    <span class="hljs-comment"># 取余，得 1</span>
<span class="hljs-built_in">print</span>(<span class="hljs-number">2</span> ** <span class="hljs-number">10</span>)  <span class="hljs-comment"># 幂运算，得 1024</span>
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">round</span>(<span class="hljs-number">3.14159</span>, <span class="hljs-number">2</span>))  <span class="hljs-comment"># 四舍五入保留 2 位</span>
</code></pre>
<h3>4.2 字符串</h3>
<pre><code class="hljs language-python">s = <span class="hljs-string">"hello"</span>
t = <span class="hljs-string">'world'</span>

<span class="hljs-comment"># 拼接与重复</span>
<span class="hljs-built_in">print</span>(s + <span class="hljs-string">" "</span> + t)   <span class="hljs-comment"># hello world</span>
<span class="hljs-built_in">print</span>(s * <span class="hljs-number">3</span>)         <span class="hljs-comment"># hellohellohello</span>

<span class="hljs-comment"># 索引与切片（左闭右开）</span>
<span class="hljs-built_in">print</span>(s[<span class="hljs-number">0</span>])          <span class="hljs-comment"># h</span>
<span class="hljs-built_in">print</span>(s[-<span class="hljs-number">1</span>])         <span class="hljs-comment"># o</span>
<span class="hljs-built_in">print</span>(s[<span class="hljs-number">1</span>:<span class="hljs-number">3</span>])        <span class="hljs-comment"># el</span>

<span class="hljs-comment"># 常用方法</span>
<span class="hljs-built_in">print</span>(s.upper())          <span class="hljs-comment"># HELLO</span>
<span class="hljs-built_in">print</span>(s.replace(<span class="hljs-string">"l"</span>, <span class="hljs-string">"L"</span>))  <span class="hljs-comment"># heLLo</span>
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">len</span>(s))             <span class="hljs-comment"># 5</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">", "</span>.join([<span class="hljs-string">"a"</span>, <span class="hljs-string">"b"</span>, <span class="hljs-string">"c"</span>]))  <span class="hljs-comment"># a, b, c</span>
</code></pre>
<h3>4.3 列表（list）</h3>
<p>列表是有序、可变的容器：</p>
<pre><code class="hljs language-python">fruits = [<span class="hljs-string">"apple"</span>, <span class="hljs-string">"banana"</span>, <span class="hljs-string">"cherry"</span>]

fruits.append(<span class="hljs-string">"orange"</span>)      <span class="hljs-comment"># 末尾追加</span>
fruits.insert(<span class="hljs-number">0</span>, <span class="hljs-string">"grape"</span>)    <span class="hljs-comment"># 指定位置插入</span>
fruits.remove(<span class="hljs-string">"apple"</span>)       <span class="hljs-comment"># 按值删除</span>
popped = fruits.pop()        <span class="hljs-comment"># 弹出末尾元素</span>

<span class="hljs-built_in">print</span>(fruits[<span class="hljs-number">0</span>])             <span class="hljs-comment"># 索引</span>
<span class="hljs-built_in">print</span>(fruits[-<span class="hljs-number">1</span>])            <span class="hljs-comment"># 最后一个</span>
<span class="hljs-built_in">print</span>(fruits[<span class="hljs-number">1</span>:<span class="hljs-number">3</span>])           <span class="hljs-comment"># 切片</span>

<span class="hljs-comment"># 列表推导式（非常常用）</span>
squares = [x ** <span class="hljs-number">2</span> <span class="hljs-keyword">for</span> x <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">10</span>)]
<span class="hljs-built_in">print</span>(squares)  <span class="hljs-comment"># [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]</span>
</code></pre>
<h3>4.4 元组（tuple）</h3>
<p>元组是不可变列表，常用于返回多个值：</p>
<pre><code class="hljs language-python">point = (<span class="hljs-number">3</span>, <span class="hljs-number">4</span>)
x, y = point          <span class="hljs-comment"># 解包</span>
<span class="hljs-built_in">print</span>(x, y)           <span class="hljs-comment"># 3 4</span>

<span class="hljs-keyword">def</span> <span class="hljs-title function_">minmax</span>(<span class="hljs-params">nums</span>):
    <span class="hljs-keyword">return</span> <span class="hljs-built_in">min</span>(nums), <span class="hljs-built_in">max</span>(nums)

lo, hi = minmax([<span class="hljs-number">3</span>, <span class="hljs-number">1</span>, <span class="hljs-number">4</span>, <span class="hljs-number">1</span>, <span class="hljs-number">5</span>])
<span class="hljs-built_in">print</span>(lo, hi)         <span class="hljs-comment"># 1 5</span>
</code></pre>
<h3>4.5 字典（dict）</h3>
<p>字典是键值对容器，查询效率高：</p>
<pre><code class="hljs language-python">person = {<span class="hljs-string">"name"</span>: <span class="hljs-string">"zorrooz"</span>, <span class="hljs-string">"age"</span>: <span class="hljs-number">25</span>}

<span class="hljs-built_in">print</span>(person[<span class="hljs-string">"name"</span>])            <span class="hljs-comment"># 取值</span>
<span class="hljs-built_in">print</span>(person.get(<span class="hljs-string">"city"</span>, <span class="hljs-string">"未知"</span>))  <span class="hljs-comment"># 安全的取值，带默认值</span>
person[<span class="hljs-string">"city"</span>] = <span class="hljs-string">"Lanzhou"</span>       <span class="hljs-comment"># 新增/修改</span>

<span class="hljs-keyword">for</span> key, value <span class="hljs-keyword">in</span> person.items():
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{key}</span>: <span class="hljs-subst">{value}</span>"</span>)

<span class="hljs-comment"># 字典推导式</span>
squares = {x: x ** <span class="hljs-number">2</span> <span class="hljs-keyword">for</span> x <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">5</span>)}
<span class="hljs-built_in">print</span>(squares)  <span class="hljs-comment"># {0: 0, 1: 1, 2: 4, 3: 9, 4: 16}</span>
</code></pre>
<h3>4.6 集合（set）</h3>
<p>集合是无序、去重的容器：</p>
<pre><code class="hljs language-python">a = {<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">3</span>, <span class="hljs-number">3</span>}
<span class="hljs-built_in">print</span>(a)              <span class="hljs-comment"># {1, 2, 3}，自动去重</span>

b = {<span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>}
<span class="hljs-built_in">print</span>(a &#x26; b)          <span class="hljs-comment"># 交集 {2, 3}</span>
<span class="hljs-built_in">print</span>(a | b)          <span class="hljs-comment"># 并集 {1, 2, 3, 4}</span>
<span class="hljs-built_in">print</span>(a - b)          <span class="hljs-comment"># 差集 {1}</span>
</code></pre>
<h2>5. 流程控制</h2>
<h3>5.1 条件判断</h3>
<pre><code class="hljs language-python">score = <span class="hljs-number">85</span>

<span class="hljs-keyword">if</span> score >= <span class="hljs-number">90</span>:
    grade = <span class="hljs-string">"A"</span>
<span class="hljs-keyword">elif</span> score >= <span class="hljs-number">80</span>:
    grade = <span class="hljs-string">"B"</span>
<span class="hljs-keyword">elif</span> score >= <span class="hljs-number">70</span>:
    grade = <span class="hljs-string">"C"</span>
<span class="hljs-keyword">else</span>:
    grade = <span class="hljs-string">"D"</span>

<span class="hljs-built_in">print</span>(<span class="hljs-string">f"等级：<span class="hljs-subst">{grade}</span>"</span>)

<span class="hljs-comment"># 三元表达式</span>
status = <span class="hljs-string">"pass"</span> <span class="hljs-keyword">if</span> score >= <span class="hljs-number">60</span> <span class="hljs-keyword">else</span> <span class="hljs-string">"fail"</span>
</code></pre>
<h3>5.2 for 循环</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 遍历可迭代对象</span>
<span class="hljs-keyword">for</span> fruit <span class="hljs-keyword">in</span> [<span class="hljs-string">"apple"</span>, <span class="hljs-string">"banana"</span>]:
    <span class="hljs-built_in">print</span>(fruit)

<span class="hljs-comment"># range 生成数列</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">5</span>):        <span class="hljs-comment"># 0 到 4</span>
    <span class="hljs-built_in">print</span>(i)

<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">1</span>, <span class="hljs-number">10</span>, <span class="hljs-number">2</span>):  <span class="hljs-comment"># 1, 3, 5, 7, 9</span>
    <span class="hljs-built_in">print</span>(i)

<span class="hljs-comment"># enumerate 同时取索引和值</span>
<span class="hljs-keyword">for</span> idx, fruit <span class="hljs-keyword">in</span> <span class="hljs-built_in">enumerate</span>([<span class="hljs-string">"a"</span>, <span class="hljs-string">"b"</span>]):
    <span class="hljs-built_in">print</span>(idx, fruit)
</code></pre>
<h3>5.3 while 循环</h3>
<pre><code class="hljs language-python">n = <span class="hljs-number">0</span>
<span class="hljs-keyword">while</span> n &#x3C; <span class="hljs-number">5</span>:
    <span class="hljs-built_in">print</span>(n)
    n += <span class="hljs-number">1</span>          <span class="hljs-comment"># 注意：忘记递增会死循环</span>
</code></pre>
<h3>5.4 break / continue / else</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># break：提前终止循环</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">10</span>):
    <span class="hljs-keyword">if</span> i == <span class="hljs-number">3</span>:
        <span class="hljs-keyword">break</span>
    <span class="hljs-built_in">print</span>(i)        <span class="hljs-comment"># 输出 0 1 2</span>

<span class="hljs-comment"># continue：跳过本次迭代</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">5</span>):
    <span class="hljs-keyword">if</span> i == <span class="hljs-number">2</span>:
        <span class="hljs-keyword">continue</span>
    <span class="hljs-built_in">print</span>(i)        <span class="hljs-comment"># 输出 0 1 3 4</span>

<span class="hljs-comment"># for-else：循环正常结束（未被 break）时执行</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">3</span>):
    <span class="hljs-keyword">if</span> i == <span class="hljs-number">99</span>:
        <span class="hljs-keyword">break</span>
<span class="hljs-keyword">else</span>:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"循环正常结束"</span>)
</code></pre>
<h2>6. 综合练习</h2>
<p>统计一段文本中每个单词的出现次数：</p>
<pre><code class="hljs language-python">text = <span class="hljs-string">"python is fun and python is powerful"</span>

words = text.split()
counter = {}

<span class="hljs-keyword">for</span> word <span class="hljs-keyword">in</span> words:
    counter[word] = counter.get(word, <span class="hljs-number">0</span>) + <span class="hljs-number">1</span>

<span class="hljs-keyword">for</span> word, count <span class="hljs-keyword">in</span> <span class="hljs-built_in">sorted</span>(counter.items(), key=<span class="hljs-keyword">lambda</span> x: -x[<span class="hljs-number">1</span>]):
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{word}</span>: <span class="hljs-subst">{count}</span>"</span>)
</code></pre>
<p>输出：</p>
<pre><code>python: 2
is: 2
fun: 1
and: 1
powerful: 1
</code></pre>
<h2>7. 小结</h2>
<ul>
<li>用 <code>venv</code> 隔离项目依赖，用 <code>pip freeze</code> 记录依赖</li>
<li>内置类型：<code>int</code> / <code>float</code> / <code>str</code> / <code>list</code> / <code>tuple</code> / <code>dict</code> / <code>set</code></li>
<li>列表与字典推导式、f-string 是日常高频语法</li>
<li>流程控制：<code>if/elif/else</code>、<code>for</code>、<code>while</code>、<code>break</code> / <code>continue</code></li>
</ul>
<p>下一篇将介绍函数、类与模块，进入真正的工程化编程。</p>`,Db=`<h1>Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas</h1>
<p>Data processing is the most core application scenario of Python. This article explains file I/O, regular expressions, and the two pillars of scientific computing, NumPy and Pandas, with examples tailored to common tasks in bioinformatics.</p>
<h2>1. File I/O</h2>
<h3>1.1 Text Files</h3>
<p>It is recommended to use the <code>with</code> statement to automatically manage file closure:</p>
<pre><code class="hljs language-python"><span class="hljs-comment"># Write</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"output.txt"</span>, <span class="hljs-string">"w"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    f.write(<span class="hljs-string">"First line\\n"</span>)
    f.write(<span class="hljs-string">"Second line\\n"</span>)

<span class="hljs-comment"># Read all</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"output.txt"</span>, <span class="hljs-string">"r"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    content = f.read()
<span class="hljs-built_in">print</span>(content)

<span class="hljs-comment"># Read line by line (recommended for large files)</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"output.txt"</span>, <span class="hljs-string">"r"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    <span class="hljs-keyword">for</span> line <span class="hljs-keyword">in</span> f:
        <span class="hljs-built_in">print</span>(line.strip())   <span class="hljs-comment"># strip removes the newline</span>
</code></pre>
<p>Mode description: <code>r</code> read, <code>w</code> write (overwrite), <code>a</code> append, <code>rb</code>/<code>wb</code> binary.</p>
<h3>1.2 Parsing FASTA Files</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">read_fasta</span>(<span class="hljs-params">path</span>):
    <span class="hljs-string">"""Read a FASTA file and return a {id: sequence} dictionary."""</span>
    sequences = {}
    current_id = <span class="hljs-literal">None</span>
    <span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(path, <span class="hljs-string">"r"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
        <span class="hljs-keyword">for</span> line <span class="hljs-keyword">in</span> f:
            line = line.strip()
            <span class="hljs-keyword">if</span> line.startswith(<span class="hljs-string">">"</span>):
                current_id = line[<span class="hljs-number">1</span>:].split()[<span class="hljs-number">0</span>]   <span class="hljs-comment"># take the first field as the id</span>
                sequences[current_id] = <span class="hljs-string">""</span>
            <span class="hljs-keyword">elif</span> current_id <span class="hljs-keyword">is</span> <span class="hljs-keyword">not</span> <span class="hljs-literal">None</span>:
                sequences[current_id] += line
    <span class="hljs-keyword">return</span> sequences

seqs = read_fasta(<span class="hljs-string">"example.fasta"</span>)
<span class="hljs-keyword">for</span> seq_id, seq <span class="hljs-keyword">in</span> seqs.items():
    <span class="hljs-built_in">print</span>(seq_id, <span class="hljs-built_in">len</span>(seq))
</code></pre>
<h3>1.3 CSV Files</h3>
<p>Standard library <code>csv</code> module:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> csv

<span class="hljs-comment"># Read</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"data.csv"</span>, <span class="hljs-string">"r"</span>, newline=<span class="hljs-string">""</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    reader = csv.DictReader(f)          <span class="hljs-comment"># read as dict by header</span>
    rows = <span class="hljs-built_in">list</span>(reader)
    <span class="hljs-built_in">print</span>(rows[<span class="hljs-number">0</span>][<span class="hljs-string">"gene"</span>], rows[<span class="hljs-number">0</span>][<span class="hljs-string">"expression"</span>])

<span class="hljs-comment"># Write</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"out.csv"</span>, <span class="hljs-string">"w"</span>, newline=<span class="hljs-string">""</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    writer = csv.writer(f)
    writer.writerow([<span class="hljs-string">"gene"</span>, <span class="hljs-string">"value"</span>])
    writer.writerows([[<span class="hljs-string">"TP53"</span>, <span class="hljs-number">12.3</span>], [<span class="hljs-string">"BRCA1"</span>, <span class="hljs-number">8.9</span>]])
</code></pre>
<blockquote>
<p>For big data scenarios, Pandas' <code>read_csv</code> is recommended (see Section 4).</p>
</blockquote>
<h2>2. Regular Expressions</h2>
<p>Regular expressions are used for text pattern matching. Common functions in the <code>re</code> module:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> re

text = <span class="hljs-string">"The protein P12345 has 500 amino acids, code: A1B2C3"</span>

<span class="hljs-comment"># match: match from the beginning</span>
m = re.<span class="hljs-keyword">match</span>(<span class="hljs-string">r"The"</span>, text)
<span class="hljs-built_in">print</span>(m.group())          <span class="hljs-comment"># The</span>

<span class="hljs-comment"># search: search for the first match</span>
m = re.search(<span class="hljs-string">r"\\d+"</span>, text)
<span class="hljs-built_in">print</span>(m.group())          <span class="hljs-comment"># 12345 (the first digit string)</span>

<span class="hljs-comment"># findall: return all matches</span>
<span class="hljs-built_in">print</span>(re.findall(<span class="hljs-string">r"[A-Z]\\d"</span>, text))   <span class="hljs-comment"># ['P1', 'A1', 'B2', 'C3']</span>

<span class="hljs-comment"># sub: substitute</span>
<span class="hljs-built_in">print</span>(re.sub(<span class="hljs-string">r"\\d+"</span>, <span class="hljs-string">"#"</span>, text))
</code></pre>
<h3>2.1 Common Metacharacters</h3>
<table>
<thead>
<tr>
<th>Pattern</th>
<th>Meaning</th>
<th>Example</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>.</code></td>
<td>Any character (except newline)</td>
<td><code>a.c</code> matches <code>abc</code></td>
</tr>
<tr>
<td><code>*</code></td>
<td>Previous character 0 or more times</td>
<td><code>ab*c</code> matches <code>ac</code>, <code>abc</code></td>
</tr>
<tr>
<td><code>+</code></td>
<td>Previous character 1 or more times</td>
<td><code>ab+c</code> matches <code>abc</code></td>
</tr>
<tr>
<td><code>?</code></td>
<td>Previous character 0 or 1 time</td>
<td><code>colou?r</code> matches color/colour</td>
</tr>
<tr>
<td><code>^</code> / <code>$</code></td>
<td>Start / end of line</td>
<td><code>^ATG</code> starts with ATG</td>
</tr>
<tr>
<td><code>[abc]</code></td>
<td>Character set</td>
<td><code>[ATGC]</code> any base</td>
</tr>
<tr>
<td><code>\\d</code> / <code>\\w</code> / <code>\\s</code></td>
<td>Digit / word character / whitespace</td>
<td></td>
</tr>
<tr>
<td><code>{n,m}</code></td>
<td>Repeat n to m times</td>
<td><code>\\d{3}</code> three digits</td>
</tr>
<tr>
<td><code>(…)</code></td>
<td>Group capture</td>
<td><code>(ATG){2}</code></td>
</tr>
</tbody>
</table>
<h3>2.2 Practical Information Extraction</h3>
<p>Extract coordinates and names from gene annotation text:</p>
<pre><code class="hljs language-python">line = <span class="hljs-string">'gene=TP53;gene_id=ENSG00000141510;chromosome=17;start=7668402;end=7687550'</span>

m = re.search(<span class="hljs-string">r"gene=(\\w+);.*?start=(\\d+);end=(\\d+)"</span>, line)
<span class="hljs-keyword">if</span> m:
    gene, start, end = m.groups()
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{gene}</span>: <span class="hljs-subst">{start}</span>-<span class="hljs-subst">{end}</span> (<span class="hljs-subst">{<span class="hljs-built_in">int</span>(end) - <span class="hljs-built_in">int</span>(start)}</span> bp)"</span>)
<span class="hljs-comment"># TP53: 7668402-7687550 (19148 bp)</span>
</code></pre>
<h3>2.3 Compiling Regular Expressions for Performance</h3>
<pre><code class="hljs language-python">motif = re.<span class="hljs-built_in">compile</span>(<span class="hljs-string">r"^ATG.*(TAA|TAG|TGA)$"</span>)   <span class="hljs-comment"># rough match for complete ORF</span>
<span class="hljs-keyword">for</span> seq <span class="hljs-keyword">in</span> [<span class="hljs-string">"ATGCCCTAA"</span>, <span class="hljs-string">"ATGGGG"</span>]:
    <span class="hljs-built_in">print</span>(seq, <span class="hljs-built_in">bool</span>(motif.<span class="hljs-keyword">match</span>(seq)))
</code></pre>
<h2>3. NumPy: Basics of Numerical Computing</h2>
<h3>3.1 Array Creation</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> numpy <span class="hljs-keyword">as</span> np

a = np.array([<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>])
b = np.zeros((<span class="hljs-number">2</span>, <span class="hljs-number">3</span>))          <span class="hljs-comment"># 2x3 all-zero matrix</span>
c = np.ones(<span class="hljs-number">5</span>)
d = np.arange(<span class="hljs-number">0</span>, <span class="hljs-number">1</span>, <span class="hljs-number">0.2</span>)      <span class="hljs-comment"># [0.  0.2 0.4 0.6 0.8]</span>
e = np.linspace(<span class="hljs-number">0</span>, <span class="hljs-number">1</span>, <span class="hljs-number">5</span>)      <span class="hljs-comment"># split 0 to 1 into 5 equal parts</span>
f = np.random.rand(<span class="hljs-number">3</span>, <span class="hljs-number">3</span>)      <span class="hljs-comment"># uniformly distributed random numbers</span>
</code></pre>
<h3>3.2 Vectorized Operations (Core Advantage)</h3>
<pre><code class="hljs language-python">x = np.array([<span class="hljs-number">1.0</span>, <span class="hljs-number">2.0</span>, <span class="hljs-number">3.0</span>])
y = np.array([<span class="hljs-number">4.0</span>, <span class="hljs-number">5.0</span>, <span class="hljs-number">6.0</span>])

<span class="hljs-built_in">print</span>(x + y)          <span class="hljs-comment"># [5. 7. 9.] without loops</span>
<span class="hljs-built_in">print</span>(x * y)          <span class="hljs-comment"># [4. 10. 18.]</span>
<span class="hljs-built_in">print</span>(np.sqrt(x))     <span class="hljs-comment"># element-wise square root</span>
<span class="hljs-built_in">print</span>(x.mean(), x.std())   <span class="hljs-comment"># 2.0 0.816...</span>

<span class="hljs-comment"># Broadcasting: automatic expansion of different shapes</span>
matrix = np.array([[<span class="hljs-number">1</span>, <span class="hljs-number">2</span>], [<span class="hljs-number">3</span>, <span class="hljs-number">4</span>]])
<span class="hljs-built_in">print</span>(matrix * <span class="hljs-number">10</span>)            <span class="hljs-comment"># multiply each element by 10</span>
<span class="hljs-built_in">print</span>(matrix + np.array([<span class="hljs-number">100</span>, <span class="hljs-number">200</span>]))  <span class="hljs-comment"># row broadcast</span>
</code></pre>
<h3>3.3 Indexing and Slicing</h3>
<pre><code class="hljs language-python">arr = np.arange(<span class="hljs-number">10</span>)
<span class="hljs-built_in">print</span>(arr[<span class="hljs-number">2</span>:<span class="hljs-number">8</span>:<span class="hljs-number">2</span>])        <span class="hljs-comment"># [2 4 6], start:stop:step</span>
<span class="hljs-built_in">print</span>(arr[arr > <span class="hljs-number">5</span>])      <span class="hljs-comment"># boolean mask [6 7 8 9]</span>

m = np.random.rand(<span class="hljs-number">4</span>, <span class="hljs-number">4</span>)
<span class="hljs-built_in">print</span>(m[<span class="hljs-number">0</span>, :])           <span class="hljs-comment"># first row</span>
<span class="hljs-built_in">print</span>(m[:, <span class="hljs-number">1</span>])           <span class="hljs-comment"># second column</span>
<span class="hljs-built_in">print</span>(m[<span class="hljs-number">1</span>:<span class="hljs-number">3</span>, <span class="hljs-number">1</span>:<span class="hljs-number">3</span>])       <span class="hljs-comment"># submatrix</span>
</code></pre>
<h3>3.4 Bioinformatics Example: Normalizing an Expression Matrix</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Simulate a gene expression matrix: 4 genes × 5 samples</span>
expr = np.array([
    [<span class="hljs-number">12.3</span>, <span class="hljs-number">15.1</span>, <span class="hljs-number">9.8</span>, <span class="hljs-number">20.4</span>, <span class="hljs-number">11.2</span>],
    [<span class="hljs-number">3.2</span>, <span class="hljs-number">4.1</span>, <span class="hljs-number">2.9</span>, <span class="hljs-number">5.0</span>, <span class="hljs-number">3.6</span>],
    [<span class="hljs-number">80.5</span>, <span class="hljs-number">92.3</span>, <span class="hljs-number">75.1</span>, <span class="hljs-number">100.2</span>, <span class="hljs-number">88.7</span>],
    [<span class="hljs-number">45.6</span>, <span class="hljs-number">41.2</span>, <span class="hljs-number">48.9</span>, <span class="hljs-number">39.8</span>, <span class="hljs-number">50.1</span>],
])

<span class="hljs-comment"># Perform Z-score normalization by row (gene): (x - mean) / std</span>
means = expr.mean(axis=<span class="hljs-number">1</span>, keepdims=<span class="hljs-literal">True</span>)
stds = expr.std(axis=<span class="hljs-number">1</span>, keepdims=<span class="hljs-literal">True</span>)
z = (expr - means) / stds

<span class="hljs-built_in">print</span>(np.<span class="hljs-built_in">round</span>(z, <span class="hljs-number">2</span>))
</code></pre>
<h2>4. Pandas: Tabular Data Processing</h2>
<h3>4.1 Series and DataFrame</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> pandas <span class="hljs-keyword">as</span> pd

<span class="hljs-comment"># Series: one-dimensional labeled array</span>
s = pd.Series([<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>], index=[<span class="hljs-string">"a"</span>, <span class="hljs-string">"b"</span>, <span class="hljs-string">"c"</span>])

<span class="hljs-comment"># DataFrame: two-dimensional table</span>
df = pd.DataFrame({
    <span class="hljs-string">"gene"</span>: [<span class="hljs-string">"TP53"</span>, <span class="hljs-string">"BRCA1"</span>, <span class="hljs-string">"EGFR"</span>, <span class="hljs-string">"MYC"</span>],
    <span class="hljs-string">"expression"</span>: [<span class="hljs-number">12.3</span>, <span class="hljs-number">8.9</span>, <span class="hljs-number">45.1</span>, <span class="hljs-number">33.2</span>],
    <span class="hljs-string">"chromosome"</span>: [<span class="hljs-number">17</span>, <span class="hljs-number">17</span>, <span class="hljs-number">7</span>, <span class="hljs-number">8</span>],
})
<span class="hljs-built_in">print</span>(df)
</code></pre>
<h3>4.2 Reading and Writing</h3>
<pre><code class="hljs language-python">df = pd.read_csv(<span class="hljs-string">"expression.csv"</span>)
df = pd.read_csv(<span class="hljs-string">"data.tsv"</span>, sep=<span class="hljs-string">"\\t"</span>)
df = pd.read_excel(<span class="hljs-string">"data.xlsx"</span>, sheet_name=<span class="hljs-string">"Sheet1"</span>)
df = pd.read_json(<span class="hljs-string">"data.json"</span>)

df.to_csv(<span class="hljs-string">"out.csv"</span>, index=<span class="hljs-literal">False</span>)
</code></pre>
<h3>4.3 Viewing and Filtering</h3>
<pre><code class="hljs language-python"><span class="hljs-built_in">print</span>(df.head())          <span class="hljs-comment"># first 5 rows</span>
<span class="hljs-built_in">print</span>(df.info())          <span class="hljs-comment"># column types and missing values</span>
<span class="hljs-built_in">print</span>(df.describe())      <span class="hljs-comment"># statistical summary of numeric columns</span>

<span class="hljs-comment"># Column selection</span>
<span class="hljs-built_in">print</span>(df[<span class="hljs-string">"gene"</span>])
<span class="hljs-built_in">print</span>(df[[<span class="hljs-string">"gene"</span>, <span class="hljs-string">"expression"</span>]])

<span class="hljs-comment"># Row filtering (boolean indexing)</span>
high = df[df[<span class="hljs-string">"expression"</span>] > <span class="hljs-number">10</span>]
<span class="hljs-built_in">print</span>(high)

<span class="hljs-comment"># Multiple conditions</span>
sel = df[(df[<span class="hljs-string">"expression"</span>] > <span class="hljs-number">10</span>) &#x26; (df[<span class="hljs-string">"chromosome"</span>] == <span class="hljs-number">17</span>)]
<span class="hljs-built_in">print</span>(sel)

<span class="hljs-comment"># isin filtering</span>
sel = df[df[<span class="hljs-string">"gene"</span>].isin([<span class="hljs-string">"TP53"</span>, <span class="hljs-string">"MYC"</span>])]
</code></pre>
<h3>4.4 Grouping and Aggregation</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Group by chromosome and compute expression mean and count</span>
grouped = df.groupby(<span class="hljs-string">"chromosome"</span>)[<span class="hljs-string">"expression"</span>].agg([<span class="hljs-string">"mean"</span>, <span class="hljs-string">"count"</span>, <span class="hljs-string">"max"</span>])
<span class="hljs-built_in">print</span>(grouped)

<span class="hljs-comment"># Transform by group (e.g., within-group standardization)</span>
df[<span class="hljs-string">"zscore"</span>] = df.groupby(<span class="hljs-string">"chromosome"</span>)[<span class="hljs-string">"expression"</span>].transform(
    <span class="hljs-keyword">lambda</span> x: (x - x.mean()) / x.std()
)
</code></pre>
<h3>4.5 Handling Missing Values</h3>
<pre><code class="hljs language-python">df = pd.DataFrame({<span class="hljs-string">"a"</span>: [<span class="hljs-number">1</span>, <span class="hljs-literal">None</span>, <span class="hljs-number">3</span>], <span class="hljs-string">"b"</span>: [<span class="hljs-number">4</span>, <span class="hljs-number">5</span>, <span class="hljs-literal">None</span>]})

<span class="hljs-built_in">print</span>(df.isna().<span class="hljs-built_in">sum</span>())        <span class="hljs-comment"># missing value count per column</span>
<span class="hljs-built_in">print</span>(df.dropna())            <span class="hljs-comment"># drop rows with missing values</span>
<span class="hljs-built_in">print</span>(df.fillna(<span class="hljs-number">0</span>))           <span class="hljs-comment"># fill with 0</span>
<span class="hljs-built_in">print</span>(df.fillna(df.mean()))   <span class="hljs-comment"># fill with mean</span>
</code></pre>
<h3>4.6 Practical Example: Filtering Differentially Expressed Genes</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Simulate: gene + treated/control expression + p-value</span>
results = pd.DataFrame({
    <span class="hljs-string">"gene"</span>: [<span class="hljs-string">f"GENE<span class="hljs-subst">{i}</span>"</span> <span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">100</span>)],
    <span class="hljs-string">"log2fc"</span>: np.random.normal(<span class="hljs-number">0</span>, <span class="hljs-number">1.5</span>, <span class="hljs-number">100</span>),
    <span class="hljs-string">"pvalue"</span>: np.random.rand(<span class="hljs-number">100</span>),
})

<span class="hljs-comment"># Filter significant genes: |log2FC| > 1 and p &#x3C; 0.05</span>
sig = results[(results[<span class="hljs-string">"log2fc"</span>].<span class="hljs-built_in">abs</span>() > <span class="hljs-number">1</span>) &#x26; (results[<span class="hljs-string">"pvalue"</span>] &#x3C; <span class="hljs-number">0.05</span>)]

<span class="hljs-comment"># Sort by log2FC</span>
sig = sig.sort_values(<span class="hljs-string">"log2fc"</span>, ascending=<span class="hljs-literal">False</span>)
<span class="hljs-built_in">print</span>(<span class="hljs-string">f"Significant gene count: <span class="hljs-subst">{<span class="hljs-built_in">len</span>(sig)}</span>"</span>)
<span class="hljs-built_in">print</span>(sig.head())
</code></pre>
<h2>5. Summary</h2>
<ul>
<li>Files: <code>with open</code> + line-by-line reading for large files; use <code>csv</code> module or Pandas for CSV</li>
<li>Regex: <code>re.search</code> / <code>findall</code> / <code>sub</code>, start with small patterns and combine</li>
<li>NumPy: vectorized operations + broadcasting, avoid Python loops</li>
<li>Pandas: <code>read_csv</code> → filtering/grouping/missing value handling → export, the standard workflow for tabular data</li>
</ul>
<p>At this point, the Python tutorial trilogy is complete: Beginner → Intermediate → Data Processing.</p>`,Ob=`<h1>Python 数据处理实战：文件 IO、正则与 NumPy/Pandas</h1>
<p>数据处理是 Python 最核心的应用场景。本文讲解文件读写、正则表达式，以及科学计算的两大支柱 NumPy 与 Pandas，示例贴合生物信息学常见任务。</p>
<h2>1. 文件读写</h2>
<h3>1.1 文本文件</h3>
<p>推荐使用 <code>with</code> 语句，自动管理文件关闭：</p>
<pre><code class="hljs language-python"><span class="hljs-comment"># 写入</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"output.txt"</span>, <span class="hljs-string">"w"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    f.write(<span class="hljs-string">"第一行\\n"</span>)
    f.write(<span class="hljs-string">"第二行\\n"</span>)

<span class="hljs-comment"># 读取全部</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"output.txt"</span>, <span class="hljs-string">"r"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    content = f.read()
<span class="hljs-built_in">print</span>(content)

<span class="hljs-comment"># 逐行读取（大文件推荐）</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"output.txt"</span>, <span class="hljs-string">"r"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    <span class="hljs-keyword">for</span> line <span class="hljs-keyword">in</span> f:
        <span class="hljs-built_in">print</span>(line.strip())   <span class="hljs-comment"># strip 去掉换行符</span>
</code></pre>
<p>模式说明：<code>r</code> 读、<code>w</code> 写（覆盖）、<code>a</code> 追加、<code>rb</code>/<code>wb</code> 二进制。</p>
<h3>1.2 解析 FASTA 文件</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">read_fasta</span>(<span class="hljs-params">path</span>):
    <span class="hljs-string">"""读取 FASTA 文件，返回 {id: 序列} 字典。"""</span>
    sequences = {}
    current_id = <span class="hljs-literal">None</span>
    <span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(path, <span class="hljs-string">"r"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
        <span class="hljs-keyword">for</span> line <span class="hljs-keyword">in</span> f:
            line = line.strip()
            <span class="hljs-keyword">if</span> line.startswith(<span class="hljs-string">">"</span>):
                current_id = line[<span class="hljs-number">1</span>:].split()[<span class="hljs-number">0</span>]   <span class="hljs-comment"># 取第一个字段作为 id</span>
                sequences[current_id] = <span class="hljs-string">""</span>
            <span class="hljs-keyword">elif</span> current_id <span class="hljs-keyword">is</span> <span class="hljs-keyword">not</span> <span class="hljs-literal">None</span>:
                sequences[current_id] += line
    <span class="hljs-keyword">return</span> sequences

seqs = read_fasta(<span class="hljs-string">"example.fasta"</span>)
<span class="hljs-keyword">for</span> seq_id, seq <span class="hljs-keyword">in</span> seqs.items():
    <span class="hljs-built_in">print</span>(seq_id, <span class="hljs-built_in">len</span>(seq))
</code></pre>
<h3>1.3 CSV 文件</h3>
<p>标准库 <code>csv</code> 模块：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> csv

<span class="hljs-comment"># 读取</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"data.csv"</span>, <span class="hljs-string">"r"</span>, newline=<span class="hljs-string">""</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    reader = csv.DictReader(f)          <span class="hljs-comment"># 按表头读取为字典</span>
    rows = <span class="hljs-built_in">list</span>(reader)
    <span class="hljs-built_in">print</span>(rows[<span class="hljs-number">0</span>][<span class="hljs-string">"gene"</span>], rows[<span class="hljs-number">0</span>][<span class="hljs-string">"expression"</span>])

<span class="hljs-comment"># 写入</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"out.csv"</span>, <span class="hljs-string">"w"</span>, newline=<span class="hljs-string">""</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    writer = csv.writer(f)
    writer.writerow([<span class="hljs-string">"gene"</span>, <span class="hljs-string">"value"</span>])
    writer.writerows([[<span class="hljs-string">"TP53"</span>, <span class="hljs-number">12.3</span>], [<span class="hljs-string">"BRCA1"</span>, <span class="hljs-number">8.9</span>]])
</code></pre>
<blockquote>
<p>大数据场景推荐 Pandas 的 <code>read_csv</code>（见第 4 节）。</p>
</blockquote>
<h2>2. 正则表达式</h2>
<p>正则用于文本模式匹配。<code>re</code> 模块常用函数：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> re

text = <span class="hljs-string">"The protein P12345 has 500 amino acids, code: A1B2C3"</span>

<span class="hljs-comment"># match：从开头匹配</span>
m = re.<span class="hljs-keyword">match</span>(<span class="hljs-string">r"The"</span>, text)
<span class="hljs-built_in">print</span>(m.group())          <span class="hljs-comment"># The</span>

<span class="hljs-comment"># search：搜索第一个匹配</span>
m = re.search(<span class="hljs-string">r"\\d+"</span>, text)
<span class="hljs-built_in">print</span>(m.group())          <span class="hljs-comment"># 12345（第一个数字串）</span>

<span class="hljs-comment"># findall：返回所有匹配</span>
<span class="hljs-built_in">print</span>(re.findall(<span class="hljs-string">r"[A-Z]\\d"</span>, text))   <span class="hljs-comment"># ['P1', 'A1', 'B2', 'C3']</span>

<span class="hljs-comment"># sub：替换</span>
<span class="hljs-built_in">print</span>(re.sub(<span class="hljs-string">r"\\d+"</span>, <span class="hljs-string">"#"</span>, text))
</code></pre>
<h3>2.1 常用元字符</h3>
<table>
<thead>
<tr>
<th>模式</th>
<th>含义</th>
<th>示例</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>.</code></td>
<td>任意字符（除换行）</td>
<td><code>a.c</code> 匹配 <code>abc</code></td>
</tr>
<tr>
<td><code>*</code></td>
<td>前一个字符 0 次或多次</td>
<td><code>ab*c</code> 匹配 <code>ac</code>、<code>abc</code></td>
</tr>
<tr>
<td><code>+</code></td>
<td>前一个字符 1 次或多次</td>
<td><code>ab+c</code> 匹配 <code>abc</code></td>
</tr>
<tr>
<td><code>?</code></td>
<td>前一个字符 0 次或 1 次</td>
<td><code>colou?r</code> 匹配 color/colour</td>
</tr>
<tr>
<td><code>^</code> / <code>$</code></td>
<td>行首 / 行尾</td>
<td><code>^ATG</code> 以 ATG 开头</td>
</tr>
<tr>
<td><code>[abc]</code></td>
<td>字符集</td>
<td><code>[ATGC]</code> 任意碱基</td>
</tr>
<tr>
<td><code>\\d</code> / <code>\\w</code> / <code>\\s</code></td>
<td>数字 / 单词字符 / 空白</td>
<td></td>
</tr>
<tr>
<td><code>{n,m}</code></td>
<td>重复 n 到 m 次</td>
<td><code>\\d{3}</code> 三个数字</td>
</tr>
<tr>
<td><code>(…)</code></td>
<td>分组捕获</td>
<td><code>(ATG){2}</code></td>
</tr>
</tbody>
</table>
<h3>2.2 提取信息实战</h3>
<p>从基因注释文本提取坐标与名称：</p>
<pre><code class="hljs language-python">line = <span class="hljs-string">'gene=TP53;gene_id=ENSG00000141510;chromosome=17;start=7668402;end=7687550'</span>

m = re.search(<span class="hljs-string">r"gene=(\\w+);.*?start=(\\d+);end=(\\d+)"</span>, line)
<span class="hljs-keyword">if</span> m:
    gene, start, end = m.groups()
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{gene}</span>: <span class="hljs-subst">{start}</span>-<span class="hljs-subst">{end}</span> (<span class="hljs-subst">{<span class="hljs-built_in">int</span>(end) - <span class="hljs-built_in">int</span>(start)}</span> bp)"</span>)
<span class="hljs-comment"># TP53: 7668402-7687550 (19148 bp)</span>
</code></pre>
<h3>2.3 编译正则提升性能</h3>
<pre><code class="hljs language-python">motif = re.<span class="hljs-built_in">compile</span>(<span class="hljs-string">r"^ATG.*(TAA|TAG|TGA)$"</span>)   <span class="hljs-comment"># 完整 ORF 粗略匹配</span>
<span class="hljs-keyword">for</span> seq <span class="hljs-keyword">in</span> [<span class="hljs-string">"ATGCCCTAA"</span>, <span class="hljs-string">"ATGGGG"</span>]:
    <span class="hljs-built_in">print</span>(seq, <span class="hljs-built_in">bool</span>(motif.<span class="hljs-keyword">match</span>(seq)))
</code></pre>
<h2>3. NumPy：数值计算基础</h2>
<h3>3.1 数组创建</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> numpy <span class="hljs-keyword">as</span> np

a = np.array([<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>])
b = np.zeros((<span class="hljs-number">2</span>, <span class="hljs-number">3</span>))          <span class="hljs-comment"># 2x3 全零矩阵</span>
c = np.ones(<span class="hljs-number">5</span>)
d = np.arange(<span class="hljs-number">0</span>, <span class="hljs-number">1</span>, <span class="hljs-number">0.2</span>)      <span class="hljs-comment"># [0.  0.2 0.4 0.6 0.8]</span>
e = np.linspace(<span class="hljs-number">0</span>, <span class="hljs-number">1</span>, <span class="hljs-number">5</span>)      <span class="hljs-comment"># 0 到 1 等分 5 份</span>
f = np.random.rand(<span class="hljs-number">3</span>, <span class="hljs-number">3</span>)      <span class="hljs-comment"># 均匀分布随机数</span>
</code></pre>
<h3>3.2 向量化运算（核心优势）</h3>
<pre><code class="hljs language-python">x = np.array([<span class="hljs-number">1.0</span>, <span class="hljs-number">2.0</span>, <span class="hljs-number">3.0</span>])
y = np.array([<span class="hljs-number">4.0</span>, <span class="hljs-number">5.0</span>, <span class="hljs-number">6.0</span>])

<span class="hljs-built_in">print</span>(x + y)          <span class="hljs-comment"># [5. 7. 9.]，无需循环</span>
<span class="hljs-built_in">print</span>(x * y)          <span class="hljs-comment"># [4. 10. 18.]</span>
<span class="hljs-built_in">print</span>(np.sqrt(x))     <span class="hljs-comment"># 逐元素开方</span>
<span class="hljs-built_in">print</span>(x.mean(), x.std())   <span class="hljs-comment"># 2.0 0.816...</span>

<span class="hljs-comment"># 广播：不同形状自动扩展</span>
matrix = np.array([[<span class="hljs-number">1</span>, <span class="hljs-number">2</span>], [<span class="hljs-number">3</span>, <span class="hljs-number">4</span>]])
<span class="hljs-built_in">print</span>(matrix * <span class="hljs-number">10</span>)            <span class="hljs-comment"># 每个元素乘 10</span>
<span class="hljs-built_in">print</span>(matrix + np.array([<span class="hljs-number">100</span>, <span class="hljs-number">200</span>]))  <span class="hljs-comment"># 行广播</span>
</code></pre>
<h3>3.3 索引与切片</h3>
<pre><code class="hljs language-python">arr = np.arange(<span class="hljs-number">10</span>)
<span class="hljs-built_in">print</span>(arr[<span class="hljs-number">2</span>:<span class="hljs-number">8</span>:<span class="hljs-number">2</span>])        <span class="hljs-comment"># [2 4 6]，起始:结束:步长</span>
<span class="hljs-built_in">print</span>(arr[arr > <span class="hljs-number">5</span>])      <span class="hljs-comment"># 布尔掩码 [6 7 8 9]</span>

m = np.random.rand(<span class="hljs-number">4</span>, <span class="hljs-number">4</span>)
<span class="hljs-built_in">print</span>(m[<span class="hljs-number">0</span>, :])           <span class="hljs-comment"># 第一行</span>
<span class="hljs-built_in">print</span>(m[:, <span class="hljs-number">1</span>])           <span class="hljs-comment"># 第二列</span>
<span class="hljs-built_in">print</span>(m[<span class="hljs-number">1</span>:<span class="hljs-number">3</span>, <span class="hljs-number">1</span>:<span class="hljs-number">3</span>])       <span class="hljs-comment"># 子矩阵</span>
</code></pre>
<h3>3.4 生物信息学示例：归一化表达矩阵</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 模拟基因表达矩阵：4 基因 × 5 样本</span>
expr = np.array([
    [<span class="hljs-number">12.3</span>, <span class="hljs-number">15.1</span>, <span class="hljs-number">9.8</span>, <span class="hljs-number">20.4</span>, <span class="hljs-number">11.2</span>],
    [<span class="hljs-number">3.2</span>, <span class="hljs-number">4.1</span>, <span class="hljs-number">2.9</span>, <span class="hljs-number">5.0</span>, <span class="hljs-number">3.6</span>],
    [<span class="hljs-number">80.5</span>, <span class="hljs-number">92.3</span>, <span class="hljs-number">75.1</span>, <span class="hljs-number">100.2</span>, <span class="hljs-number">88.7</span>],
    [<span class="hljs-number">45.6</span>, <span class="hljs-number">41.2</span>, <span class="hljs-number">48.9</span>, <span class="hljs-number">39.8</span>, <span class="hljs-number">50.1</span>],
])

<span class="hljs-comment"># 按行（基因）做 Z-score 标准化：(x - mean) / std</span>
means = expr.mean(axis=<span class="hljs-number">1</span>, keepdims=<span class="hljs-literal">True</span>)
stds = expr.std(axis=<span class="hljs-number">1</span>, keepdims=<span class="hljs-literal">True</span>)
z = (expr - means) / stds

<span class="hljs-built_in">print</span>(np.<span class="hljs-built_in">round</span>(z, <span class="hljs-number">2</span>))
</code></pre>
<h2>4. Pandas：表格数据处理</h2>
<h3>4.1 Series 与 DataFrame</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> pandas <span class="hljs-keyword">as</span> pd

<span class="hljs-comment"># Series：一维带标签数组</span>
s = pd.Series([<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>], index=[<span class="hljs-string">"a"</span>, <span class="hljs-string">"b"</span>, <span class="hljs-string">"c"</span>])

<span class="hljs-comment"># DataFrame：二维表格</span>
df = pd.DataFrame({
    <span class="hljs-string">"gene"</span>: [<span class="hljs-string">"TP53"</span>, <span class="hljs-string">"BRCA1"</span>, <span class="hljs-string">"EGFR"</span>, <span class="hljs-string">"MYC"</span>],
    <span class="hljs-string">"expression"</span>: [<span class="hljs-number">12.3</span>, <span class="hljs-number">8.9</span>, <span class="hljs-number">45.1</span>, <span class="hljs-number">33.2</span>],
    <span class="hljs-string">"chromosome"</span>: [<span class="hljs-number">17</span>, <span class="hljs-number">17</span>, <span class="hljs-number">7</span>, <span class="hljs-number">8</span>],
})
<span class="hljs-built_in">print</span>(df)
</code></pre>
<h3>4.2 读取与写入</h3>
<pre><code class="hljs language-python">df = pd.read_csv(<span class="hljs-string">"expression.csv"</span>)
df = pd.read_csv(<span class="hljs-string">"data.tsv"</span>, sep=<span class="hljs-string">"\\t"</span>)
df = pd.read_excel(<span class="hljs-string">"data.xlsx"</span>, sheet_name=<span class="hljs-string">"Sheet1"</span>)
df = pd.read_json(<span class="hljs-string">"data.json"</span>)

df.to_csv(<span class="hljs-string">"out.csv"</span>, index=<span class="hljs-literal">False</span>)
</code></pre>
<h3>4.3 查看与筛选</h3>
<pre><code class="hljs language-python"><span class="hljs-built_in">print</span>(df.head())          <span class="hljs-comment"># 前 5 行</span>
<span class="hljs-built_in">print</span>(df.info())          <span class="hljs-comment"># 列类型与缺失值</span>
<span class="hljs-built_in">print</span>(df.describe())      <span class="hljs-comment"># 数值列统计摘要</span>

<span class="hljs-comment"># 列选择</span>
<span class="hljs-built_in">print</span>(df[<span class="hljs-string">"gene"</span>])
<span class="hljs-built_in">print</span>(df[[<span class="hljs-string">"gene"</span>, <span class="hljs-string">"expression"</span>]])

<span class="hljs-comment"># 行筛选（布尔索引）</span>
high = df[df[<span class="hljs-string">"expression"</span>] > <span class="hljs-number">10</span>]
<span class="hljs-built_in">print</span>(high)

<span class="hljs-comment"># 多条件</span>
sel = df[(df[<span class="hljs-string">"expression"</span>] > <span class="hljs-number">10</span>) &#x26; (df[<span class="hljs-string">"chromosome"</span>] == <span class="hljs-number">17</span>)]
<span class="hljs-built_in">print</span>(sel)

<span class="hljs-comment"># isin 筛选</span>
sel = df[df[<span class="hljs-string">"gene"</span>].isin([<span class="hljs-string">"TP53"</span>, <span class="hljs-string">"MYC"</span>])]
</code></pre>
<h3>4.4 分组与聚合</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 按染色体分组，统计表达均值与数量</span>
grouped = df.groupby(<span class="hljs-string">"chromosome"</span>)[<span class="hljs-string">"expression"</span>].agg([<span class="hljs-string">"mean"</span>, <span class="hljs-string">"count"</span>, <span class="hljs-string">"max"</span>])
<span class="hljs-built_in">print</span>(grouped)

<span class="hljs-comment"># 按分组转换（如组内标准化）</span>
df[<span class="hljs-string">"zscore"</span>] = df.groupby(<span class="hljs-string">"chromosome"</span>)[<span class="hljs-string">"expression"</span>].transform(
    <span class="hljs-keyword">lambda</span> x: (x - x.mean()) / x.std()
)
</code></pre>
<h3>4.5 缺失值处理</h3>
<pre><code class="hljs language-python">df = pd.DataFrame({<span class="hljs-string">"a"</span>: [<span class="hljs-number">1</span>, <span class="hljs-literal">None</span>, <span class="hljs-number">3</span>], <span class="hljs-string">"b"</span>: [<span class="hljs-number">4</span>, <span class="hljs-number">5</span>, <span class="hljs-literal">None</span>]})

<span class="hljs-built_in">print</span>(df.isna().<span class="hljs-built_in">sum</span>())        <span class="hljs-comment"># 每列缺失数量</span>
<span class="hljs-built_in">print</span>(df.dropna())            <span class="hljs-comment"># 删除含缺失的行</span>
<span class="hljs-built_in">print</span>(df.fillna(<span class="hljs-number">0</span>))           <span class="hljs-comment"># 填充为 0</span>
<span class="hljs-built_in">print</span>(df.fillna(df.mean()))   <span class="hljs-comment"># 用均值填充</span>
</code></pre>
<h3>4.6 实战：差异表达基因筛选</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 模拟：基因 + 处理组/对照组表达 + p 值</span>
results = pd.DataFrame({
    <span class="hljs-string">"gene"</span>: [<span class="hljs-string">f"GENE<span class="hljs-subst">{i}</span>"</span> <span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">100</span>)],
    <span class="hljs-string">"log2fc"</span>: np.random.normal(<span class="hljs-number">0</span>, <span class="hljs-number">1.5</span>, <span class="hljs-number">100</span>),
    <span class="hljs-string">"pvalue"</span>: np.random.rand(<span class="hljs-number">100</span>),
})

<span class="hljs-comment"># 筛选显著差异基因：|log2FC| > 1 且 p &#x3C; 0.05</span>
sig = results[(results[<span class="hljs-string">"log2fc"</span>].<span class="hljs-built_in">abs</span>() > <span class="hljs-number">1</span>) &#x26; (results[<span class="hljs-string">"pvalue"</span>] &#x3C; <span class="hljs-number">0.05</span>)]

<span class="hljs-comment"># 按 log2FC 排序</span>
sig = sig.sort_values(<span class="hljs-string">"log2fc"</span>, ascending=<span class="hljs-literal">False</span>)
<span class="hljs-built_in">print</span>(<span class="hljs-string">f"显著基因数: <span class="hljs-subst">{<span class="hljs-built_in">len</span>(sig)}</span>"</span>)
<span class="hljs-built_in">print</span>(sig.head())
</code></pre>
<h2>5. 小结</h2>
<ul>
<li>文件：<code>with open</code> + 逐行读取处理大文件；CSV 用 <code>csv</code> 模块或 Pandas</li>
<li>正则：<code>re.search</code> / <code>findall</code> / <code>sub</code>，先写小模式再组合</li>
<li>NumPy：向量化运算 + 广播，避免 Python 循环</li>
<li>Pandas：<code>read_csv</code> → 筛选/分组/缺失值处理 → 导出，是表格数据的标准工作流</li>
</ul>
<p>至此 Python 教程三部曲完成：入门 → 进阶 → 数据处理。</p>`,Ib=`<h1>R Language Primer: Data Structures, Vectorization, and Functions</h1>
<p>R is the go-to language for statistical computing and data visualization. This article explains the core of R: five data structures, vectorized operations, control flow, and functions, laying the groundwork for subsequent tidyverse and ggplot2.</p>
<h2>1. Environment Setup</h2>
<h3>1.1 Installation and IDE</h3>
<ul>
<li>Install <a href="https://cran.r-project.org/">R</a></li>
<li>Recommended IDE: <a href="https://posit.co/download/rstudio-desktop/">RStudio</a> (free), or VS Code + R extension</li>
</ul>
<h3>1.2 Basic Operations</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Assignment: &#x3C;- is the classic R assignment operator (= also works)</span>
x <span class="hljs-operator">&#x3C;-</span> 10
x

<span class="hljs-comment"># View help</span>
<span class="hljs-operator">?</span>mean
help<span class="hljs-punctuation">(</span><span class="hljs-string">"lm"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>2. Basic Data Structures</h2>
<p>R's core data structures can be classified by dimension and homogeneity:</p>
<table>
<thead>
<tr>
<th>Structure</th>
<th>Dimension</th>
<th>Element Type</th>
</tr>
</thead>
<tbody>
<tr>
<td>vector</td>
<td>1</td>
<td>Homogeneous</td>
</tr>
<tr>
<td>factor</td>
<td>1</td>
<td>Categorical</td>
</tr>
<tr>
<td>matrix</td>
<td>2</td>
<td>Homogeneous</td>
</tr>
<tr>
<td>data.frame</td>
<td>2</td>
<td>Heterogeneous allowed</td>
</tr>
<tr>
<td>list</td>
<td>N</td>
<td>Heterogeneous allowed</td>
</tr>
</tbody>
</table>
<h3>2.1 Vectors</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Create vectors</span>
a <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span>            <span class="hljs-comment"># combine</span>
b <span class="hljs-operator">&#x3C;-</span> 1<span class="hljs-operator">:</span><span class="hljs-number">10</span>                        <span class="hljs-comment"># sequence</span>
<span class="hljs-built_in">c</span> <span class="hljs-operator">&#x3C;-</span> seq<span class="hljs-punctuation">(</span><span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-number">0.2</span><span class="hljs-punctuation">)</span>         <span class="hljs-comment"># step sequence</span>
d <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-string">"A"</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span>                 <span class="hljs-comment"># repeat</span>
e <span class="hljs-operator">&#x3C;-</span> sample<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">100</span><span class="hljs-punctuation">,</span> <span class="hljs-number">10</span><span class="hljs-punctuation">)</span>           <span class="hljs-comment"># random sample</span>

<span class="hljs-comment"># Types</span>
typeof<span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span>      <span class="hljs-comment"># "double"</span>
<span class="hljs-built_in">is.numeric</span><span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span>  <span class="hljs-comment"># TRUE</span>

<span class="hljs-comment"># Indexing (R starts at 1!)</span>
a<span class="hljs-punctuation">[</span><span class="hljs-number">1</span><span class="hljs-punctuation">]</span>           <span class="hljs-comment"># 1</span>
a<span class="hljs-punctuation">[</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">]</span>     <span class="hljs-comment"># elements 1 and 3</span>
a<span class="hljs-punctuation">[</span><span class="hljs-operator">-</span><span class="hljs-number">2</span><span class="hljs-punctuation">]</span>          <span class="hljs-comment"># remove element 2</span>
a<span class="hljs-punctuation">[</span>a <span class="hljs-operator">></span> <span class="hljs-number">3</span><span class="hljs-punctuation">]</span>       <span class="hljs-comment"># conditional filtering</span>
<span class="hljs-built_in">names</span><span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span> <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"x1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x2"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x3"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x4"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x5"</span><span class="hljs-punctuation">)</span>
a<span class="hljs-punctuation">[</span><span class="hljs-string">"x1"</span><span class="hljs-punctuation">]</span>        <span class="hljs-comment"># access by name</span>
</code></pre>
<blockquote>
<p><strong>R indexes from 1</strong>, which is the biggest habit difference from Python.</p>
</blockquote>
<h3>2.2 Factors</h3>
<pre><code class="hljs language-r">treatment <span class="hljs-operator">&#x3C;-</span> factor<span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"drug"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"drug"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
                    levels <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"drug"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
levels<span class="hljs-punctuation">(</span>treatment<span class="hljs-punctuation">)</span>          <span class="hljs-comment"># "control" "drug"</span>
table<span class="hljs-punctuation">(</span>treatment<span class="hljs-punctuation">)</span>           <span class="hljs-comment"># frequency counts</span>
</code></pre>
<h3>2.3 Matrices</h3>
<pre><code class="hljs language-r">m <span class="hljs-operator">&#x3C;-</span> matrix<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">12</span><span class="hljs-punctuation">,</span> nrow <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> ncol <span class="hljs-operator">=</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span>
<span class="hljs-comment">#      [,1] [,2] [,3] [,4]</span>
<span class="hljs-comment"># [1,]    1    4    7   10</span>
<span class="hljs-comment"># [2,]    2    5    8   11</span>
<span class="hljs-comment"># [3,]    3    6    9   12</span>

m<span class="hljs-punctuation">[</span><span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">]</span>        <span class="hljs-comment"># row 2, column 3: 8</span>
m<span class="hljs-punctuation">[</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-punctuation">]</span>         <span class="hljs-comment"># row 1</span>
m<span class="hljs-punctuation">[</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">]</span>         <span class="hljs-comment"># column 2</span>
t<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">)</span>           <span class="hljs-comment"># transpose</span>
m <span class="hljs-operator">%*%</span> t<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">)</span>     <span class="hljs-comment"># matrix multiplication</span>
</code></pre>
<h3>2.4 Data Frames</h3>
<p>Data frames are the core structure for tabular data:</p>
<pre><code class="hljs language-r">df <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  <span class="hljs-built_in">expression</span> <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">12.3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.9</span><span class="hljs-punctuation">,</span> <span class="hljs-number">45.1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  chromosome <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># View structure</span>
str<span class="hljs-punctuation">(</span>df<span class="hljs-punctuation">)</span>
summary<span class="hljs-punctuation">(</span>df<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Extract columns</span>
df<span class="hljs-operator">$</span><span class="hljs-built_in">expression</span>          <span class="hljs-comment"># $ operator</span>
df<span class="hljs-punctuation">[[</span><span class="hljs-string">"expression"</span><span class="hljs-punctuation">]</span><span class="hljs-punctuation">]</span>
df<span class="hljs-punctuation">[</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"expression"</span><span class="hljs-punctuation">]</span>

<span class="hljs-comment"># Extract rows</span>
df<span class="hljs-punctuation">[</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-punctuation">]</span>                <span class="hljs-comment"># first row</span>
df<span class="hljs-punctuation">[</span>df<span class="hljs-operator">$</span><span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">10</span><span class="hljs-punctuation">,</span> <span class="hljs-punctuation">]</span>

<span class="hljs-comment"># Add column</span>
df<span class="hljs-operator">$</span>log2fc <span class="hljs-operator">&#x3C;-</span> log2<span class="hljs-punctuation">(</span>df<span class="hljs-operator">$</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.5 Lists</h3>
<p>Lists can hold elements of any type:</p>
<pre><code class="hljs language-r">result <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span>
  name <span class="hljs-operator">=</span> <span class="hljs-string">"Analysis Result"</span><span class="hljs-punctuation">,</span>
  pvalue <span class="hljs-operator">=</span> <span class="hljs-number">0.003</span><span class="hljs-punctuation">,</span>
  coef <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1.2</span><span class="hljs-punctuation">,</span> <span class="hljs-operator">-</span><span class="hljs-number">0.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3.1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  model <span class="hljs-operator">=</span> lm<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span> <span class="hljs-operator">~</span> chromosome<span class="hljs-punctuation">,</span> data <span class="hljs-operator">=</span> df<span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

result<span class="hljs-operator">$</span>pvalue
result<span class="hljs-punctuation">[[</span><span class="hljs-number">3</span><span class="hljs-punctuation">]</span><span class="hljs-punctuation">]</span>            <span class="hljs-comment"># numeric vector</span>
</code></pre>
<h2>3. Vectorized Operations</h2>
<p><strong>R's core idea: operate on entire vectors, not element-by-element loops.</strong></p>
<pre><code class="hljs language-r">x <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span>

x <span class="hljs-operator">*</span> <span class="hljs-number">2</span>              <span class="hljs-comment"># element-wise multiplication</span>
x <span class="hljs-operator">+</span> <span class="hljs-number">10</span>
<span class="hljs-built_in">sqrt</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
log2<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># vector with vector</span>
y <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">)</span>
x <span class="hljs-operator">+</span> y              <span class="hljs-comment"># add corresponding positions</span>

<span class="hljs-comment"># comparison operations return logical vectors</span>
x <span class="hljs-operator">></span> <span class="hljs-number">3</span>              <span class="hljs-comment"># FALSE FALSE FALSE  TRUE  TRUE</span>

<span class="hljs-comment"># statistical functions</span>
<span class="hljs-built_in">sum</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; mean<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; median<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; sd<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; <span class="hljs-built_in">range</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
</code></pre>
<h2>4. Control Flow</h2>
<h3>4.1 Conditionals</h3>
<pre><code class="hljs language-r">score <span class="hljs-operator">&#x3C;-</span> 85

<span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>score <span class="hljs-operator">>=</span> <span class="hljs-number">90</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  grade <span class="hljs-operator">&#x3C;-</span> <span class="hljs-string">"A"</span>
<span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>score <span class="hljs-operator">>=</span> <span class="hljs-number">80</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  grade <span class="hljs-operator">&#x3C;-</span> <span class="hljs-string">"B"</span>
<span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-punctuation">{</span>
  grade <span class="hljs-operator">&#x3C;-</span> <span class="hljs-string">"C"</span>
<span class="hljs-punctuation">}</span>
print<span class="hljs-punctuation">(</span>grade<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># vectorized conditional: ifelse</span>
values <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">60</span><span class="hljs-punctuation">,</span> <span class="hljs-number">90</span><span class="hljs-punctuation">,</span> <span class="hljs-number">45</span><span class="hljs-punctuation">,</span> <span class="hljs-number">75</span><span class="hljs-punctuation">)</span>
result <span class="hljs-operator">&#x3C;-</span> ifelse<span class="hljs-punctuation">(</span>values <span class="hljs-operator">>=</span> <span class="hljs-number">60</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"pass"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"fail"</span><span class="hljs-punctuation">)</span>
print<span class="hljs-punctuation">(</span>result<span class="hljs-punctuation">)</span>   <span class="hljs-comment"># "pass" "pass" "fail" "pass"</span>
</code></pre>
<h3>4.2 Loops</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># for loop</span>
<span class="hljs-keyword">for</span> <span class="hljs-punctuation">(</span>i <span class="hljs-keyword">in</span> <span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">5</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  print<span class="hljs-punctuation">(</span>i <span class="hljs-operator">^</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">}</span>

<span class="hljs-comment"># iterate over vector elements</span>
genes <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">)</span>
<span class="hljs-keyword">for</span> <span class="hljs-punctuation">(</span>g <span class="hljs-keyword">in</span> genes<span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  print<span class="hljs-punctuation">(</span>paste<span class="hljs-punctuation">(</span><span class="hljs-string">"Analyzing"</span><span class="hljs-punctuation">,</span> g<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">}</span>

<span class="hljs-comment"># while loop</span>
n <span class="hljs-operator">&#x3C;-</span> 0
<span class="hljs-keyword">while</span> <span class="hljs-punctuation">(</span>n <span class="hljs-operator">&#x3C;</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  n <span class="hljs-operator">&#x3C;-</span> n <span class="hljs-operator">+</span> <span class="hljs-number">1</span>
<span class="hljs-punctuation">}</span>
</code></pre>
<h3>4.3 Avoid Loops When Possible: The apply Family</h3>
<p>R commonly uses the <code>apply</code> family instead of explicit loops:</p>
<pre><code class="hljs language-r">m <span class="hljs-operator">&#x3C;-</span> matrix<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">12</span><span class="hljs-punctuation">,</span> nrow <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span>

apply<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> mean<span class="hljs-punctuation">)</span>      <span class="hljs-comment"># row-wise means</span>
apply<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-built_in">sum</span><span class="hljs-punctuation">)</span>       <span class="hljs-comment"># column-wise sums</span>

lapply<span class="hljs-punctuation">(</span><span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-operator">:</span><span class="hljs-number">6</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> mean<span class="hljs-punctuation">)</span>   <span class="hljs-comment"># apply function to each element of a list, returns a list</span>
sapply<span class="hljs-punctuation">(</span><span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-operator">:</span><span class="hljs-number">6</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> mean<span class="hljs-punctuation">)</span>   <span class="hljs-comment"># simplified version, returns a vector</span>
</code></pre>
<h2>5. Functions</h2>
<h3>5.1 Defining Functions</h3>
<pre><code class="hljs language-r">gc_content <span class="hljs-operator">&#x3C;-</span> <span class="hljs-keyword">function</span><span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  seq <span class="hljs-operator">&#x3C;-</span> toupper<span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">)</span>
  gc <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">sum</span><span class="hljs-punctuation">(</span>strsplit<span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">,</span> <span class="hljs-string">""</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">[[</span><span class="hljs-number">1</span><span class="hljs-punctuation">]</span><span class="hljs-punctuation">]</span> <span class="hljs-operator">%in%</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"G"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"C"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
  gc <span class="hljs-operator">/</span> nchar<span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">)</span> <span class="hljs-operator">*</span> <span class="hljs-number">100</span>
<span class="hljs-punctuation">}</span>

gc_content<span class="hljs-punctuation">(</span><span class="hljs-string">"ATGCCGA"</span><span class="hljs-punctuation">)</span>
<span class="hljs-comment"># [1] 57.14286</span>
</code></pre>
<h3>5.2 Default Arguments and Return Values</h3>
<pre><code class="hljs language-r">normalize <span class="hljs-operator">&#x3C;-</span> <span class="hljs-keyword">function</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">,</span> method <span class="hljs-operator">=</span> <span class="hljs-string">"zscore"</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  <span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>method <span class="hljs-operator">==</span> <span class="hljs-string">"zscore"</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
    <span class="hljs-punctuation">(</span>x <span class="hljs-operator">-</span> mean<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">/</span> sd<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>method <span class="hljs-operator">==</span> <span class="hljs-string">"minmax"</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
    <span class="hljs-punctuation">(</span>x <span class="hljs-operator">-</span> <span class="hljs-built_in">min</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">/</span> <span class="hljs-punctuation">(</span><span class="hljs-built_in">max</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span> <span class="hljs-operator">-</span> <span class="hljs-built_in">min</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-punctuation">{</span>
    stop<span class="hljs-punctuation">(</span><span class="hljs-string">"Unknown method: "</span><span class="hljs-punctuation">,</span> method<span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">}</span>
<span class="hljs-punctuation">}</span>

normalize<span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
normalize<span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> method <span class="hljs-operator">=</span> <span class="hljs-string">"minmax"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>5.3 Anonymous Functions and Functional Programming</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># anonymous functions</span>
sapply<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">5</span><span class="hljs-punctuation">,</span> <span class="hljs-keyword">function</span><span class="hljs-punctuation">(</span>i<span class="hljs-punctuation">)</span> i <span class="hljs-operator">^</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># purrr-style map (requires tidyverse, detailed in the next article)</span>
<span class="hljs-comment"># purrr::map_dbl(1:5, ~ .x ^ 2)</span>
</code></pre>
<h2>6. Practical Tips</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># pipe (native pipe |>/ in R 4.1+ or magrittr %>%)</span>
x <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span>
x <span class="hljs-operator">|></span> mean<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> <span class="hljs-built_in">round</span><span class="hljs-punctuation">(</span><span class="hljs-number">2</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># handle missing values</span>
x <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-literal">NA</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-literal">NA</span><span class="hljs-punctuation">)</span>
<span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
<span class="hljs-built_in">sum</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
x<span class="hljs-punctuation">[</span><span class="hljs-operator">!</span><span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">]</span>
na.omit<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># merging data frames</span>
df1 <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>id <span class="hljs-operator">=</span> <span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> value <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"a"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"b"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"c"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
df2 <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>id <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-operator">:</span><span class="hljs-number">4</span><span class="hljs-punctuation">,</span> score <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">90</span><span class="hljs-punctuation">,</span> <span class="hljs-number">85</span><span class="hljs-punctuation">,</span> <span class="hljs-number">88</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
merged <span class="hljs-operator">&#x3C;-</span> merge<span class="hljs-punctuation">(</span>df1<span class="hljs-punctuation">,</span> df2<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"id"</span><span class="hljs-punctuation">)</span>   <span class="hljs-comment"># similar to SQL join</span>
</code></pre>
<h2>7. Summary</h2>
<ul>
<li>Five core structures: vectors / factors / matrices / data frames / lists</li>
<li><strong>R indexing starts at 1</strong>; <code>$</code> for columns; <code>[ ]</code> for rows</li>
<li>Vectorized operations + <code>ifelse</code> + the <code>apply</code> family replace explicit loops</li>
<li>Functions: default arguments, <code>stop()</code> for errors, anonymous functions</li>
</ul>
<p>The next article will introduce tidyverse: using <code>dplyr</code> for elegant data manipulation.</p>`,Nb=`<h1>R 语言入门：数据结构、向量化与函数</h1>
<p>R 是统计计算与数据可视化的首选语言。本文讲解 R 的核心：五种数据结构、向量化运算思维、控制流与函数，为后续 tidyverse 与 ggplot2 打下基础。</p>
<h2>1. 环境准备</h2>
<h3>1.1 安装与 IDE</h3>
<ul>
<li>安装 <a href="https://cran.r-project.org/">R</a></li>
<li>推荐 IDE：<a href="https://posit.co/download/rstudio-desktop/">RStudio</a>（免费），或 VS Code + R 扩展</li>
</ul>
<h3>1.2 基本操作</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 赋值：&#x3C;- 是 R 的经典赋值符（= 也可用）</span>
x <span class="hljs-operator">&#x3C;-</span> 10
x

<span class="hljs-comment"># 查看帮助</span>
<span class="hljs-operator">?</span>mean
help<span class="hljs-punctuation">(</span><span class="hljs-string">"lm"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>2. 基础数据结构</h2>
<p>R 的核心数据结构按维度与同质性划分：</p>
<table>
<thead>
<tr>
<th>结构</th>
<th>维度</th>
<th>元素类型</th>
</tr>
</thead>
<tbody>
<tr>
<td>vector（向量）</td>
<td>1</td>
<td>同质</td>
</tr>
<tr>
<td>factor（因子）</td>
<td>1</td>
<td>分类</td>
</tr>
<tr>
<td>matrix（矩阵）</td>
<td>2</td>
<td>同质</td>
</tr>
<tr>
<td>data.frame（数据框）</td>
<td>2</td>
<td>可异质</td>
</tr>
<tr>
<td>list（列表）</td>
<td>N</td>
<td>可异质</td>
</tr>
</tbody>
</table>
<h3>2.1 向量（vector）</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 创建向量</span>
a <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span>            <span class="hljs-comment"># combine</span>
b <span class="hljs-operator">&#x3C;-</span> 1<span class="hljs-operator">:</span><span class="hljs-number">10</span>                        <span class="hljs-comment"># 序列</span>
<span class="hljs-built_in">c</span> <span class="hljs-operator">&#x3C;-</span> seq<span class="hljs-punctuation">(</span><span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-number">0.2</span><span class="hljs-punctuation">)</span>         <span class="hljs-comment"># 步长序列</span>
d <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-string">"A"</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span>                 <span class="hljs-comment"># 重复</span>
e <span class="hljs-operator">&#x3C;-</span> sample<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">100</span><span class="hljs-punctuation">,</span> <span class="hljs-number">10</span><span class="hljs-punctuation">)</span>           <span class="hljs-comment"># 随机抽样</span>

<span class="hljs-comment"># 类型</span>
typeof<span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span>      <span class="hljs-comment"># "double"</span>
<span class="hljs-built_in">is.numeric</span><span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span>  <span class="hljs-comment"># TRUE</span>

<span class="hljs-comment"># 索引（R 从 1 开始！）</span>
a<span class="hljs-punctuation">[</span><span class="hljs-number">1</span><span class="hljs-punctuation">]</span>           <span class="hljs-comment"># 1</span>
a<span class="hljs-punctuation">[</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">]</span>     <span class="hljs-comment"># 第 1、3 个元素</span>
a<span class="hljs-punctuation">[</span><span class="hljs-operator">-</span><span class="hljs-number">2</span><span class="hljs-punctuation">]</span>          <span class="hljs-comment"># 删除第 2 个</span>
a<span class="hljs-punctuation">[</span>a <span class="hljs-operator">></span> <span class="hljs-number">3</span><span class="hljs-punctuation">]</span>       <span class="hljs-comment"># 条件筛选</span>
<span class="hljs-built_in">names</span><span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span> <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"x1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x2"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x3"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x4"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x5"</span><span class="hljs-punctuation">)</span>
a<span class="hljs-punctuation">[</span><span class="hljs-string">"x1"</span><span class="hljs-punctuation">]</span>        <span class="hljs-comment"># 按名字取</span>
</code></pre>
<blockquote>
<p><strong>R 的索引从 1 开始</strong>，这是与 Python 最大的习惯差异。</p>
</blockquote>
<h3>2.2 因子（factor）</h3>
<pre><code class="hljs language-r">treatment <span class="hljs-operator">&#x3C;-</span> factor<span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"drug"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"drug"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
                    levels <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"drug"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
levels<span class="hljs-punctuation">(</span>treatment<span class="hljs-punctuation">)</span>          <span class="hljs-comment"># "control" "drug"</span>
table<span class="hljs-punctuation">(</span>treatment<span class="hljs-punctuation">)</span>           <span class="hljs-comment"># 频数统计</span>
</code></pre>
<h3>2.3 矩阵（matrix）</h3>
<pre><code class="hljs language-r">m <span class="hljs-operator">&#x3C;-</span> matrix<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">12</span><span class="hljs-punctuation">,</span> nrow <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> ncol <span class="hljs-operator">=</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span>
<span class="hljs-comment">#      [,1] [,2] [,3] [,4]</span>
<span class="hljs-comment"># [1,]    1    4    7   10</span>
<span class="hljs-comment"># [2,]    2    5    8   11</span>
<span class="hljs-comment"># [3,]    3    6    9   12</span>

m<span class="hljs-punctuation">[</span><span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">]</span>        <span class="hljs-comment"># 第 2 行第 3 列：8</span>
m<span class="hljs-punctuation">[</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-punctuation">]</span>         <span class="hljs-comment"># 第 1 行</span>
m<span class="hljs-punctuation">[</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">]</span>         <span class="hljs-comment"># 第 2 列</span>
t<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">)</span>           <span class="hljs-comment"># 转置</span>
m <span class="hljs-operator">%*%</span> t<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">)</span>     <span class="hljs-comment"># 矩阵乘法</span>
</code></pre>
<h3>2.4 数据框（data.frame）</h3>
<p>数据框是表格数据的核心结构：</p>
<pre><code class="hljs language-r">df <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  <span class="hljs-built_in">expression</span> <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">12.3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.9</span><span class="hljs-punctuation">,</span> <span class="hljs-number">45.1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  chromosome <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 查看结构</span>
str<span class="hljs-punctuation">(</span>df<span class="hljs-punctuation">)</span>
summary<span class="hljs-punctuation">(</span>df<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 取列</span>
df<span class="hljs-operator">$</span><span class="hljs-built_in">expression</span>          <span class="hljs-comment"># $ 运算符</span>
df<span class="hljs-punctuation">[[</span><span class="hljs-string">"expression"</span><span class="hljs-punctuation">]</span><span class="hljs-punctuation">]</span>
df<span class="hljs-punctuation">[</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"expression"</span><span class="hljs-punctuation">]</span>

<span class="hljs-comment"># 取行</span>
df<span class="hljs-punctuation">[</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-punctuation">]</span>                <span class="hljs-comment"># 第一行</span>
df<span class="hljs-punctuation">[</span>df<span class="hljs-operator">$</span><span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">10</span><span class="hljs-punctuation">,</span> <span class="hljs-punctuation">]</span>

<span class="hljs-comment"># 增加列</span>
df<span class="hljs-operator">$</span>log2fc <span class="hljs-operator">&#x3C;-</span> log2<span class="hljs-punctuation">(</span>df<span class="hljs-operator">$</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.5 列表（list）</h3>
<p>列表可容纳任意类型：</p>
<pre><code class="hljs language-r">result <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span>
  name <span class="hljs-operator">=</span> <span class="hljs-string">"分析结果"</span><span class="hljs-punctuation">,</span>
  pvalue <span class="hljs-operator">=</span> <span class="hljs-number">0.003</span><span class="hljs-punctuation">,</span>
  coef <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1.2</span><span class="hljs-punctuation">,</span> <span class="hljs-operator">-</span><span class="hljs-number">0.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3.1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  model <span class="hljs-operator">=</span> lm<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span> <span class="hljs-operator">~</span> chromosome<span class="hljs-punctuation">,</span> data <span class="hljs-operator">=</span> df<span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

result<span class="hljs-operator">$</span>pvalue
result<span class="hljs-punctuation">[[</span><span class="hljs-number">3</span><span class="hljs-punctuation">]</span><span class="hljs-punctuation">]</span>            <span class="hljs-comment"># 数值向量</span>
</code></pre>
<h2>3. 向量化运算</h2>
<p><strong>R 的核心思想：对整个向量运算，而非逐元素循环。</strong></p>
<pre><code class="hljs language-r">x <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span>

x <span class="hljs-operator">*</span> <span class="hljs-number">2</span>              <span class="hljs-comment"># 逐元素乘</span>
x <span class="hljs-operator">+</span> <span class="hljs-number">10</span>
<span class="hljs-built_in">sqrt</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
log2<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 向量与向量</span>
y <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">)</span>
x <span class="hljs-operator">+</span> y              <span class="hljs-comment"># 对应位置相加</span>

<span class="hljs-comment"># 比较运算返回逻辑向量</span>
x <span class="hljs-operator">></span> <span class="hljs-number">3</span>              <span class="hljs-comment"># FALSE FALSE FALSE  TRUE  TRUE</span>

<span class="hljs-comment"># 统计函数</span>
<span class="hljs-built_in">sum</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; mean<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; median<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; sd<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; <span class="hljs-built_in">range</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
</code></pre>
<h2>4. 控制流</h2>
<h3>4.1 条件</h3>
<pre><code class="hljs language-r">score <span class="hljs-operator">&#x3C;-</span> 85

<span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>score <span class="hljs-operator">>=</span> <span class="hljs-number">90</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  grade <span class="hljs-operator">&#x3C;-</span> <span class="hljs-string">"A"</span>
<span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>score <span class="hljs-operator">>=</span> <span class="hljs-number">80</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  grade <span class="hljs-operator">&#x3C;-</span> <span class="hljs-string">"B"</span>
<span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-punctuation">{</span>
  grade <span class="hljs-operator">&#x3C;-</span> <span class="hljs-string">"C"</span>
<span class="hljs-punctuation">}</span>
print<span class="hljs-punctuation">(</span>grade<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 向量化条件：ifelse</span>
values <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">60</span><span class="hljs-punctuation">,</span> <span class="hljs-number">90</span><span class="hljs-punctuation">,</span> <span class="hljs-number">45</span><span class="hljs-punctuation">,</span> <span class="hljs-number">75</span><span class="hljs-punctuation">)</span>
result <span class="hljs-operator">&#x3C;-</span> ifelse<span class="hljs-punctuation">(</span>values <span class="hljs-operator">>=</span> <span class="hljs-number">60</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"pass"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"fail"</span><span class="hljs-punctuation">)</span>
print<span class="hljs-punctuation">(</span>result<span class="hljs-punctuation">)</span>   <span class="hljs-comment"># "pass" "pass" "fail" "pass"</span>
</code></pre>
<h3>4.2 循环</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># for 循环</span>
<span class="hljs-keyword">for</span> <span class="hljs-punctuation">(</span>i <span class="hljs-keyword">in</span> <span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">5</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  print<span class="hljs-punctuation">(</span>i <span class="hljs-operator">^</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">}</span>

<span class="hljs-comment"># 遍历向量元素</span>
genes <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">)</span>
<span class="hljs-keyword">for</span> <span class="hljs-punctuation">(</span>g <span class="hljs-keyword">in</span> genes<span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  print<span class="hljs-punctuation">(</span>paste<span class="hljs-punctuation">(</span><span class="hljs-string">"分析"</span><span class="hljs-punctuation">,</span> g<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">}</span>

<span class="hljs-comment"># while 循环</span>
n <span class="hljs-operator">&#x3C;-</span> 0
<span class="hljs-keyword">while</span> <span class="hljs-punctuation">(</span>n <span class="hljs-operator">&#x3C;</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  n <span class="hljs-operator">&#x3C;-</span> n <span class="hljs-operator">+</span> <span class="hljs-number">1</span>
<span class="hljs-punctuation">}</span>
</code></pre>
<h3>4.3 尽量避免循环：apply 家族</h3>
<p>R 中常用 <code>apply</code> 家族替代显式循环：</p>
<pre><code class="hljs language-r">m <span class="hljs-operator">&#x3C;-</span> matrix<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">12</span><span class="hljs-punctuation">,</span> nrow <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span>

apply<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> mean<span class="hljs-punctuation">)</span>      <span class="hljs-comment"># 对每行求均值</span>
apply<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-built_in">sum</span><span class="hljs-punctuation">)</span>       <span class="hljs-comment"># 对每列求和</span>

lapply<span class="hljs-punctuation">(</span><span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-operator">:</span><span class="hljs-number">6</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> mean<span class="hljs-punctuation">)</span>   <span class="hljs-comment"># 对列表每个元素操作，返回列表</span>
sapply<span class="hljs-punctuation">(</span><span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-operator">:</span><span class="hljs-number">6</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> mean<span class="hljs-punctuation">)</span>   <span class="hljs-comment"># 简化版，返回向量</span>
</code></pre>
<h2>5. 函数</h2>
<h3>5.1 定义函数</h3>
<pre><code class="hljs language-r">gc_content <span class="hljs-operator">&#x3C;-</span> <span class="hljs-keyword">function</span><span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  seq <span class="hljs-operator">&#x3C;-</span> toupper<span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">)</span>
  gc <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">sum</span><span class="hljs-punctuation">(</span>strsplit<span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">,</span> <span class="hljs-string">""</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">[[</span><span class="hljs-number">1</span><span class="hljs-punctuation">]</span><span class="hljs-punctuation">]</span> <span class="hljs-operator">%in%</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"G"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"C"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
  gc <span class="hljs-operator">/</span> nchar<span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">)</span> <span class="hljs-operator">*</span> <span class="hljs-number">100</span>
<span class="hljs-punctuation">}</span>

gc_content<span class="hljs-punctuation">(</span><span class="hljs-string">"ATGCCGA"</span><span class="hljs-punctuation">)</span>
<span class="hljs-comment"># [1] 57.14286</span>
</code></pre>
<h3>5.2 默认参数与返回值</h3>
<pre><code class="hljs language-r">normalize <span class="hljs-operator">&#x3C;-</span> <span class="hljs-keyword">function</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">,</span> method <span class="hljs-operator">=</span> <span class="hljs-string">"zscore"</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  <span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>method <span class="hljs-operator">==</span> <span class="hljs-string">"zscore"</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
    <span class="hljs-punctuation">(</span>x <span class="hljs-operator">-</span> mean<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">/</span> sd<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>method <span class="hljs-operator">==</span> <span class="hljs-string">"minmax"</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
    <span class="hljs-punctuation">(</span>x <span class="hljs-operator">-</span> <span class="hljs-built_in">min</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">/</span> <span class="hljs-punctuation">(</span><span class="hljs-built_in">max</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span> <span class="hljs-operator">-</span> <span class="hljs-built_in">min</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-punctuation">{</span>
    stop<span class="hljs-punctuation">(</span><span class="hljs-string">"未知方法: "</span><span class="hljs-punctuation">,</span> method<span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">}</span>
<span class="hljs-punctuation">}</span>

normalize<span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
normalize<span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> method <span class="hljs-operator">=</span> <span class="hljs-string">"minmax"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>5.3 匿名函数与函数式编程</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 匿名函数</span>
sapply<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">5</span><span class="hljs-punctuation">,</span> <span class="hljs-keyword">function</span><span class="hljs-punctuation">(</span>i<span class="hljs-punctuation">)</span> i <span class="hljs-operator">^</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># purrr 风格的 map（需要 tidyverse，下一篇详述）</span>
<span class="hljs-comment"># purrr::map_dbl(1:5, ~ .x ^ 2)</span>
</code></pre>
<h2>6. 实用技巧</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># 管道（R 4.1+ 原生管道 |>/ 或 magrittr %>%）</span>
x <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span>
x <span class="hljs-operator">|></span> mean<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> <span class="hljs-built_in">round</span><span class="hljs-punctuation">(</span><span class="hljs-number">2</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 缺失值处理</span>
x <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-literal">NA</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-literal">NA</span><span class="hljs-punctuation">)</span>
<span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
<span class="hljs-built_in">sum</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
x<span class="hljs-punctuation">[</span><span class="hljs-operator">!</span><span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">]</span>
na.omit<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 合并数据框</span>
df1 <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>id <span class="hljs-operator">=</span> <span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> value <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"a"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"b"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"c"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
df2 <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>id <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-operator">:</span><span class="hljs-number">4</span><span class="hljs-punctuation">,</span> score <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">90</span><span class="hljs-punctuation">,</span> <span class="hljs-number">85</span><span class="hljs-punctuation">,</span> <span class="hljs-number">88</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
merged <span class="hljs-operator">&#x3C;-</span> merge<span class="hljs-punctuation">(</span>df1<span class="hljs-punctuation">,</span> df2<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"id"</span><span class="hljs-punctuation">)</span>   <span class="hljs-comment"># 类似 SQL join</span>
</code></pre>
<h2>7. 小结</h2>
<ul>
<li>五种核心结构：向量 / 因子 / 矩阵 / 数据框 / 列表</li>
<li><strong>R 索引从 1 开始</strong>；<code>$</code> 取列；<code>[ ]</code> 取行</li>
<li>向量化运算 + <code>ifelse</code> + <code>apply</code> 家族替代显式循环</li>
<li>函数：默认参数、<code>stop()</code> 报错、匿名函数</li>
</ul>
<p>下一篇将介绍 tidyverse：用 <code>dplyr</code> 优雅地完成数据操作。</p>`,Fb=`<h1>R ggplot2 Data Visualization: Grammar of Graphics and Common Charts</h1>
<p>ggplot2 is R's most powerful plotting package, built on the <strong>Grammar of Graphics</strong>: data, mapping, geometric objects, statistical transformations, coordinates, facets, and themes, stacked layer by layer into a figure.</p>
<h2>1. Basic Principles</h2>
<pre><code class="hljs language-r">library<span class="hljs-punctuation">(</span>ggplot2<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Skeleton of a plot</span>
ggplot<span class="hljs-punctuation">(</span>data <span class="hljs-operator">=</span> dataset<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> variable<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> variable<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_xxx<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>          <span class="hljs-comment"># Geometric layer (point/line/bar/boxplot...)</span>
  labs<span class="hljs-punctuation">(</span>...<span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>           <span class="hljs-comment"># Labels</span>
  theme_xxx<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>           <span class="hljs-comment"># Theme</span>
</code></pre>
<ul>
<li><code>data</code>: data frame</li>
<li><code>aes()</code>: aesthetic mapping (x, y, color, fill, size, shape)</li>
<li><code>geom_*</code>: geometric objects, one per layer</li>
</ul>
<h2>2. Quick Start with Built-in Datasets</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># mpg: car fuel economy data; iris: iris flower data</span>
head<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">)</span>
head<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.1 Scatter Plot</h3>
<pre><code class="hljs language-r">ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Color by category + smooth trend line</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-built_in">class</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> alpha <span class="hljs-operator">=</span> <span class="hljs-number">0.7</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_smooth<span class="hljs-punctuation">(</span>method <span class="hljs-operator">=</span> <span class="hljs-string">"lm"</span><span class="hljs-punctuation">,</span> se <span class="hljs-operator">=</span> <span class="hljs-literal">FALSE</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"black"</span><span class="hljs-punctuation">,</span> linewidth <span class="hljs-operator">=</span> <span class="hljs-number">0.6</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.2 Boxplot</h3>
<pre><code class="hljs language-r">ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> fill <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_boxplot<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_jitter<span class="hljs-punctuation">(</span>width <span class="hljs-operator">=</span> <span class="hljs-number">0.2</span><span class="hljs-punctuation">,</span> size <span class="hljs-operator">=</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> alpha <span class="hljs-operator">=</span> <span class="hljs-number">0.5</span><span class="hljs-punctuation">)</span>   <span class="hljs-comment"># Overlay points to show distribution</span>
</code></pre>
<h3>2.3 Histogram</h3>
<pre><code class="hljs language-r">ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Petal.Length<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_histogram<span class="hljs-punctuation">(</span>bins <span class="hljs-operator">=</span> <span class="hljs-number">30</span><span class="hljs-punctuation">,</span> fill <span class="hljs-operator">=</span> <span class="hljs-string">"steelblue"</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"white"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.4 Bar Chart</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Count bar chart</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> <span class="hljs-built_in">class</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_bar<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Value bar chart (requires aggregation first)</span>
library<span class="hljs-punctuation">(</span>dplyr<span class="hljs-punctuation">)</span>
avg_hwy <span class="hljs-operator">&#x3C;-</span> mpg <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span><span class="hljs-built_in">class</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>mean_hwy <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span>hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span><span class="hljs-punctuation">)</span>

ggplot<span class="hljs-punctuation">(</span>avg_hwy<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> reorder<span class="hljs-punctuation">(</span><span class="hljs-built_in">class</span><span class="hljs-punctuation">,</span> mean_hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> mean_hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_col<span class="hljs-punctuation">(</span>fill <span class="hljs-operator">=</span> <span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  coord_flip<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>          <span class="hljs-comment"># Horizontal, easier to read long category names</span>
</code></pre>
<h3>2.5 Line Chart</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Time series example</span>
df <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>
  time <span class="hljs-operator">=</span> seq<span class="hljs-punctuation">(</span><span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">20</span><span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  value <span class="hljs-operator">=</span> <span class="hljs-built_in">sin</span><span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">(</span><span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">20</span><span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span> rnorm<span class="hljs-punctuation">(</span><span class="hljs-number">11</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0.1</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

ggplot<span class="hljs-punctuation">(</span>df<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> time<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> value<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_line<span class="hljs-punctuation">(</span>linewidth <span class="hljs-operator">=</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>3. Faceting (facet)</h2>
<p>Split into multiple subplots by variables:</p>
<pre><code class="hljs language-r"><span class="hljs-comment"># Facet by drv: one row</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  facet_wrap<span class="hljs-punctuation">(</span><span class="hljs-operator">~</span>drv<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Two-variable faceting: grid</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  facet_grid<span class="hljs-punctuation">(</span>cyl <span class="hljs-operator">~</span> drv<span class="hljs-punctuation">)</span>
</code></pre>
<h2>4. Color Themes</h2>
<h3>4.1 Scale Control</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Continuous color scale</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> year<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  scale_color_gradient<span class="hljs-punctuation">(</span>low <span class="hljs-operator">=</span> <span class="hljs-string">"blue"</span><span class="hljs-punctuation">,</span> high <span class="hljs-operator">=</span> <span class="hljs-string">"red"</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Discrete color palette</span>
ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Width<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2.5</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  scale_color_brewer<span class="hljs-punctuation">(</span>palette <span class="hljs-operator">=</span> <span class="hljs-string">"Set1"</span><span class="hljs-punctuation">)</span>     <span class="hljs-comment"># Built-in RColorBrewer palette</span>

<span class="hljs-comment"># Manually specified</span>
scale_color_manual<span class="hljs-punctuation">(</span>values <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"#10b981"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"#f59e0b"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.2 Theme Customization</h3>
<pre><code class="hljs language-r">p <span class="hljs-operator">&#x3C;-</span> ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Width<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2.5</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Classic white background theme</span>
p <span class="hljs-operator">+</span> theme_bw<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Minimal theme</span>
p <span class="hljs-operator">+</span> theme_minimal<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Custom details</span>
p <span class="hljs-operator">+</span>
  theme_minimal<span class="hljs-punctuation">(</span>base_size <span class="hljs-operator">=</span> <span class="hljs-number">14</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  theme<span class="hljs-punctuation">(</span>
    legend.position <span class="hljs-operator">=</span> <span class="hljs-string">"top"</span><span class="hljs-punctuation">,</span>                 <span class="hljs-comment"># Legend position</span>
    panel.grid.minor <span class="hljs-operator">=</span> element_blank<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>      <span class="hljs-comment"># Remove minor grid lines</span>
    plot.title <span class="hljs-operator">=</span> element_text<span class="hljs-punctuation">(</span>face <span class="hljs-operator">=</span> <span class="hljs-string">"bold"</span><span class="hljs-punctuation">,</span> size <span class="hljs-operator">=</span> <span class="hljs-number">16</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">)</span>
</code></pre>
<h2>5. Labels and Export</h2>
<pre><code class="hljs language-r">p <span class="hljs-operator">&#x3C;-</span> ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Width<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2.5</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  labs<span class="hljs-punctuation">(</span>
    title <span class="hljs-operator">=</span> <span class="hljs-string">"Iris Measurement Data"</span><span class="hljs-punctuation">,</span>
    subtitle <span class="hljs-operator">=</span> <span class="hljs-string">"Sepal size distribution"</span><span class="hljs-punctuation">,</span>
    x <span class="hljs-operator">=</span> <span class="hljs-string">"Sepal Length (cm)"</span><span class="hljs-punctuation">,</span>
    y <span class="hljs-operator">=</span> <span class="hljs-string">"Sepal Width (cm)"</span><span class="hljs-punctuation">,</span>
    color <span class="hljs-operator">=</span> <span class="hljs-string">"Species"</span>
  <span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  theme_minimal<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Save: ggsave automatically recognizes the extension</span>
ggsave<span class="hljs-punctuation">(</span><span class="hljs-string">"iris_scatter.png"</span><span class="hljs-punctuation">,</span> p<span class="hljs-punctuation">,</span> width <span class="hljs-operator">=</span> <span class="hljs-number">8</span><span class="hljs-punctuation">,</span> height <span class="hljs-operator">=</span> <span class="hljs-number">6</span><span class="hljs-punctuation">,</span> dpi <span class="hljs-operator">=</span> <span class="hljs-number">300</span><span class="hljs-punctuation">)</span>
ggsave<span class="hljs-punctuation">(</span><span class="hljs-string">"iris_scatter.pdf"</span><span class="hljs-punctuation">,</span> p<span class="hljs-punctuation">,</span> width <span class="hljs-operator">=</span> <span class="hljs-number">8</span><span class="hljs-punctuation">,</span> height <span class="hljs-operator">=</span> <span class="hljs-number">6</span><span class="hljs-punctuation">)</span>
</code></pre>
<blockquote>
<p>For journal submission, export <strong>PDF vector graphics</strong> (300 dpi PNG for web).</p>
</blockquote>
<h2>6. Bioinformatics Practice: Volcano Plot</h2>
<p>Standard visualization for differential expression analysis:</p>
<pre><code class="hljs language-r"><span class="hljs-comment"># Simulate differential expression results</span>
set.seed<span class="hljs-punctuation">(</span><span class="hljs-number">42</span><span class="hljs-punctuation">)</span>
de <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> paste0<span class="hljs-punctuation">(</span><span class="hljs-string">"GENE"</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">500</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  log2fc <span class="hljs-operator">=</span> rnorm<span class="hljs-punctuation">(</span><span class="hljs-number">500</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1.8</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  pvalue <span class="hljs-operator">=</span> runif<span class="hljs-punctuation">(</span><span class="hljs-number">500</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Mark significant genes</span>
de <span class="hljs-operator">&#x3C;-</span> de <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>significance <span class="hljs-operator">=</span> case_when<span class="hljs-punctuation">(</span>
    <span class="hljs-built_in">abs</span><span class="hljs-punctuation">(</span>log2fc<span class="hljs-punctuation">)</span> <span class="hljs-operator">></span> <span class="hljs-number">1</span> <span class="hljs-operator">&#x26;</span> pvalue <span class="hljs-operator">&#x3C;</span> <span class="hljs-number">0.05</span> <span class="hljs-operator">~</span> <span class="hljs-string">"up/down"</span><span class="hljs-punctuation">,</span>
    pvalue <span class="hljs-operator">&#x3C;</span> <span class="hljs-number">0.05</span> <span class="hljs-operator">~</span> <span class="hljs-string">"significant"</span><span class="hljs-punctuation">,</span>
    <span class="hljs-literal">TRUE</span> <span class="hljs-operator">~</span> <span class="hljs-string">"ns"</span>
  <span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>

ggplot<span class="hljs-punctuation">(</span>de<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> log2fc<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> <span class="hljs-operator">-</span>log10<span class="hljs-punctuation">(</span>pvalue<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> significance<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">1.5</span><span class="hljs-punctuation">,</span> alpha <span class="hljs-operator">=</span> <span class="hljs-number">0.6</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  scale_color_manual<span class="hljs-punctuation">(</span>values <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"up/down"</span> <span class="hljs-operator">=</span> <span class="hljs-string">"#ef4444"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"significant"</span> <span class="hljs-operator">=</span> <span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"ns"</span> <span class="hljs-operator">=</span> <span class="hljs-string">"#9ca3af"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_vline<span class="hljs-punctuation">(</span>xintercept <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-operator">-</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> linetype <span class="hljs-operator">=</span> <span class="hljs-string">"dashed"</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"gray50"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_hline<span class="hljs-punctuation">(</span>yintercept <span class="hljs-operator">=</span> <span class="hljs-operator">-</span>log10<span class="hljs-punctuation">(</span><span class="hljs-number">0.05</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> linetype <span class="hljs-operator">=</span> <span class="hljs-string">"dashed"</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"gray50"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  labs<span class="hljs-punctuation">(</span>title <span class="hljs-operator">=</span> <span class="hljs-string">"Differential Expression Volcano Plot"</span><span class="hljs-punctuation">,</span> x <span class="hljs-operator">=</span> <span class="hljs-string">"log2 Fold Change"</span><span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> <span class="hljs-string">"-log10(p-value)"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  theme_minimal<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>7. Summary</h2>
<ul>
<li>Grammar of Graphics: <code>ggplot(data, aes()) + geom_xxx() + labs() + theme_xxx()</code></li>
<li>Common geometries: <code>geom_point</code> / <code>geom_boxplot</code> / <code>geom_histogram</code> / <code>geom_col</code> / <code>geom_line</code></li>
<li>Group by <code>aes(color/fill)</code>, split plots by <code>facet_wrap</code> / <code>facet_grid</code></li>
<li><code>theme_minimal</code> + <code>labs</code> + <code>ggsave(dpi=300)</code> delivers publication-ready output</li>
</ul>
<p>With this, the R tutorial trilogy is complete: Introduction → tidyverse → ggplot2.</p>`,$b=`<h1>R ggplot2 数据可视化：图层语法与常用图表</h1>
<p>ggplot2 是 R 最强大的绘图包，基于<strong>图层语法</strong>（Grammar of Graphics）：数据、映射、几何对象、统计变换、坐标、分面、主题，逐层叠加成图。</p>
<h2>1. 基本原理</h2>
<pre><code class="hljs language-r">library<span class="hljs-punctuation">(</span>ggplot2<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 一张图的骨架</span>
ggplot<span class="hljs-punctuation">(</span>data <span class="hljs-operator">=</span> 数据集<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> 变量<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> 变量<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_xxx<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>          <span class="hljs-comment"># 几何图层（点/线/柱/箱线...）</span>
  labs<span class="hljs-punctuation">(</span>...<span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>           <span class="hljs-comment"># 标签</span>
  theme_xxx<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>           <span class="hljs-comment"># 主题</span>
</code></pre>
<ul>
<li><code>data</code>：数据框</li>
<li><code>aes()</code>：美学映射（x、y、color、fill、size、shape）</li>
<li><code>geom_*</code>：几何对象，每层一个</li>
</ul>
<h2>2. 内置数据集快速上手</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># mpg：汽车油耗数据；iris：鸢尾花数据</span>
head<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">)</span>
head<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.1 散点图（scatter）</h3>
<pre><code class="hljs language-r">ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 按类别着色 + 平滑趋势线</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-built_in">class</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> alpha <span class="hljs-operator">=</span> <span class="hljs-number">0.7</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_smooth<span class="hljs-punctuation">(</span>method <span class="hljs-operator">=</span> <span class="hljs-string">"lm"</span><span class="hljs-punctuation">,</span> se <span class="hljs-operator">=</span> <span class="hljs-literal">FALSE</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"black"</span><span class="hljs-punctuation">,</span> linewidth <span class="hljs-operator">=</span> <span class="hljs-number">0.6</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.2 箱线图（boxplot）</h3>
<pre><code class="hljs language-r">ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> fill <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_boxplot<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_jitter<span class="hljs-punctuation">(</span>width <span class="hljs-operator">=</span> <span class="hljs-number">0.2</span><span class="hljs-punctuation">,</span> size <span class="hljs-operator">=</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> alpha <span class="hljs-operator">=</span> <span class="hljs-number">0.5</span><span class="hljs-punctuation">)</span>   <span class="hljs-comment"># 叠加散点展示分布</span>
</code></pre>
<h3>2.3 直方图（histogram）</h3>
<pre><code class="hljs language-r">ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Petal.Length<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_histogram<span class="hljs-punctuation">(</span>bins <span class="hljs-operator">=</span> <span class="hljs-number">30</span><span class="hljs-punctuation">,</span> fill <span class="hljs-operator">=</span> <span class="hljs-string">"steelblue"</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"white"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.4 柱状图（bar）</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 计数柱状图</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> <span class="hljs-built_in">class</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_bar<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 数值柱状图（需要先聚合）</span>
library<span class="hljs-punctuation">(</span>dplyr<span class="hljs-punctuation">)</span>
avg_hwy <span class="hljs-operator">&#x3C;-</span> mpg <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span><span class="hljs-built_in">class</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>mean_hwy <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span>hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span><span class="hljs-punctuation">)</span>

ggplot<span class="hljs-punctuation">(</span>avg_hwy<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> reorder<span class="hljs-punctuation">(</span><span class="hljs-built_in">class</span><span class="hljs-punctuation">,</span> mean_hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> mean_hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_col<span class="hljs-punctuation">(</span>fill <span class="hljs-operator">=</span> <span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  coord_flip<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>          <span class="hljs-comment"># 横置，长类别名更易读</span>
</code></pre>
<h3>2.5 折线图（line）</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 时间序列示例</span>
df <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>
  time <span class="hljs-operator">=</span> seq<span class="hljs-punctuation">(</span><span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">20</span><span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  value <span class="hljs-operator">=</span> <span class="hljs-built_in">sin</span><span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">(</span><span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">20</span><span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span> rnorm<span class="hljs-punctuation">(</span><span class="hljs-number">11</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0.1</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

ggplot<span class="hljs-punctuation">(</span>df<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> time<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> value<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_line<span class="hljs-punctuation">(</span>linewidth <span class="hljs-operator">=</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>3. 分面（facet）</h2>
<p>按变量拆分成多个子图：</p>
<pre><code class="hljs language-r"><span class="hljs-comment"># 按 drv 分面：一排</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  facet_wrap<span class="hljs-punctuation">(</span><span class="hljs-operator">~</span>drv<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 双变量分面：网格</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  facet_grid<span class="hljs-punctuation">(</span>cyl <span class="hljs-operator">~</span> drv<span class="hljs-punctuation">)</span>
</code></pre>
<h2>4. 颜色主题</h2>
<h3>4.1 标度控制（scale）</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 连续色阶</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> year<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  scale_color_gradient<span class="hljs-punctuation">(</span>low <span class="hljs-operator">=</span> <span class="hljs-string">"blue"</span><span class="hljs-punctuation">,</span> high <span class="hljs-operator">=</span> <span class="hljs-string">"red"</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 离散色板</span>
ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Width<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2.5</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  scale_color_brewer<span class="hljs-punctuation">(</span>palette <span class="hljs-operator">=</span> <span class="hljs-string">"Set1"</span><span class="hljs-punctuation">)</span>     <span class="hljs-comment"># 内置 RColorBrewer 色板</span>

<span class="hljs-comment"># 手动指定</span>
scale_color_manual<span class="hljs-punctuation">(</span>values <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"#10b981"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"#f59e0b"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.2 主题美化</h3>
<pre><code class="hljs language-r">p <span class="hljs-operator">&#x3C;-</span> ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Width<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2.5</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 经典白底主题</span>
p <span class="hljs-operator">+</span> theme_bw<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 极简主题</span>
p <span class="hljs-operator">+</span> theme_minimal<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 自定义细节</span>
p <span class="hljs-operator">+</span>
  theme_minimal<span class="hljs-punctuation">(</span>base_size <span class="hljs-operator">=</span> <span class="hljs-number">14</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  theme<span class="hljs-punctuation">(</span>
    legend.position <span class="hljs-operator">=</span> <span class="hljs-string">"top"</span><span class="hljs-punctuation">,</span>                 <span class="hljs-comment"># 图例位置</span>
    panel.grid.minor <span class="hljs-operator">=</span> element_blank<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>      <span class="hljs-comment"># 去掉次要网格线</span>
    plot.title <span class="hljs-operator">=</span> element_text<span class="hljs-punctuation">(</span>face <span class="hljs-operator">=</span> <span class="hljs-string">"bold"</span><span class="hljs-punctuation">,</span> size <span class="hljs-operator">=</span> <span class="hljs-number">16</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">)</span>
</code></pre>
<h2>5. 标签与导出</h2>
<pre><code class="hljs language-r">p <span class="hljs-operator">&#x3C;-</span> ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Width<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2.5</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  labs<span class="hljs-punctuation">(</span>
    title <span class="hljs-operator">=</span> <span class="hljs-string">"鸢尾花测量数据"</span><span class="hljs-punctuation">,</span>
    subtitle <span class="hljs-operator">=</span> <span class="hljs-string">"Sepal 尺寸分布"</span><span class="hljs-punctuation">,</span>
    x <span class="hljs-operator">=</span> <span class="hljs-string">"花萼长度 (cm)"</span><span class="hljs-punctuation">,</span>
    y <span class="hljs-operator">=</span> <span class="hljs-string">"花萼宽度 (cm)"</span><span class="hljs-punctuation">,</span>
    color <span class="hljs-operator">=</span> <span class="hljs-string">"物种"</span>
  <span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  theme_minimal<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 保存：ggsave 自动识别扩展名</span>
ggsave<span class="hljs-punctuation">(</span><span class="hljs-string">"iris_scatter.png"</span><span class="hljs-punctuation">,</span> p<span class="hljs-punctuation">,</span> width <span class="hljs-operator">=</span> <span class="hljs-number">8</span><span class="hljs-punctuation">,</span> height <span class="hljs-operator">=</span> <span class="hljs-number">6</span><span class="hljs-punctuation">,</span> dpi <span class="hljs-operator">=</span> <span class="hljs-number">300</span><span class="hljs-punctuation">)</span>
ggsave<span class="hljs-punctuation">(</span><span class="hljs-string">"iris_scatter.pdf"</span><span class="hljs-punctuation">,</span> p<span class="hljs-punctuation">,</span> width <span class="hljs-operator">=</span> <span class="hljs-number">8</span><span class="hljs-punctuation">,</span> height <span class="hljs-operator">=</span> <span class="hljs-number">6</span><span class="hljs-punctuation">)</span>
</code></pre>
<blockquote>
<p>论文投稿建议导出 <strong>PDF 矢量图</strong>（300 dpi 的 PNG 用于网页）。</p>
</blockquote>
<h2>6. 生物信息学实战：火山图</h2>
<p>差异表达分析的标准可视化：</p>
<pre><code class="hljs language-r"><span class="hljs-comment"># 模拟差异表达结果</span>
set.seed<span class="hljs-punctuation">(</span><span class="hljs-number">42</span><span class="hljs-punctuation">)</span>
de <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> paste0<span class="hljs-punctuation">(</span><span class="hljs-string">"GENE"</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">500</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  log2fc <span class="hljs-operator">=</span> rnorm<span class="hljs-punctuation">(</span><span class="hljs-number">500</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1.8</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  pvalue <span class="hljs-operator">=</span> runif<span class="hljs-punctuation">(</span><span class="hljs-number">500</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 标记显著基因</span>
de <span class="hljs-operator">&#x3C;-</span> de <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>significance <span class="hljs-operator">=</span> case_when<span class="hljs-punctuation">(</span>
    <span class="hljs-built_in">abs</span><span class="hljs-punctuation">(</span>log2fc<span class="hljs-punctuation">)</span> <span class="hljs-operator">></span> <span class="hljs-number">1</span> <span class="hljs-operator">&#x26;</span> pvalue <span class="hljs-operator">&#x3C;</span> <span class="hljs-number">0.05</span> <span class="hljs-operator">~</span> <span class="hljs-string">"up/down"</span><span class="hljs-punctuation">,</span>
    pvalue <span class="hljs-operator">&#x3C;</span> <span class="hljs-number">0.05</span> <span class="hljs-operator">~</span> <span class="hljs-string">"significant"</span><span class="hljs-punctuation">,</span>
    <span class="hljs-literal">TRUE</span> <span class="hljs-operator">~</span> <span class="hljs-string">"ns"</span>
  <span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>

ggplot<span class="hljs-punctuation">(</span>de<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> log2fc<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> <span class="hljs-operator">-</span>log10<span class="hljs-punctuation">(</span>pvalue<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> significance<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">1.5</span><span class="hljs-punctuation">,</span> alpha <span class="hljs-operator">=</span> <span class="hljs-number">0.6</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  scale_color_manual<span class="hljs-punctuation">(</span>values <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"up/down"</span> <span class="hljs-operator">=</span> <span class="hljs-string">"#ef4444"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"significant"</span> <span class="hljs-operator">=</span> <span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"ns"</span> <span class="hljs-operator">=</span> <span class="hljs-string">"#9ca3af"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_vline<span class="hljs-punctuation">(</span>xintercept <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-operator">-</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> linetype <span class="hljs-operator">=</span> <span class="hljs-string">"dashed"</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"gray50"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_hline<span class="hljs-punctuation">(</span>yintercept <span class="hljs-operator">=</span> <span class="hljs-operator">-</span>log10<span class="hljs-punctuation">(</span><span class="hljs-number">0.05</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> linetype <span class="hljs-operator">=</span> <span class="hljs-string">"dashed"</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"gray50"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  labs<span class="hljs-punctuation">(</span>title <span class="hljs-operator">=</span> <span class="hljs-string">"差异表达火山图"</span><span class="hljs-punctuation">,</span> x <span class="hljs-operator">=</span> <span class="hljs-string">"log2 Fold Change"</span><span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> <span class="hljs-string">"-log10(p-value)"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  theme_minimal<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>7. 小结</h2>
<ul>
<li>图层语法：<code>ggplot(data, aes()) + geom_xxx() + labs() + theme_xxx()</code></li>
<li>常用几何：<code>geom_point</code> / <code>geom_boxplot</code> / <code>geom_histogram</code> / <code>geom_col</code> / <code>geom_line</code></li>
<li>分组用 <code>aes(color/fill)</code>，拆图用 <code>facet_wrap</code> / <code>facet_grid</code></li>
<li>主题 <code>theme_minimal</code> + <code>labs</code> + <code>ggsave(dpi=300)</code> 即出版级输出</li>
</ul>
<p>至此 R 教程三部曲完成：入门 → tidyverse → ggplot2。</p>`,Bb=`<h1>R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning</h1>
<p>tidyverse is a collection of R packages with a consistent design philosophy, where <code>dplyr</code> handles data manipulation and <code>tidyr</code> handles data reshaping. The core idea is <strong>chaining verbs with pipes</strong>, making every step of data processing clear and readable.</p>
<h2>1. Installation and Loading</h2>
<pre><code class="hljs language-r">install.packages<span class="hljs-punctuation">(</span><span class="hljs-string">"tidyverse"</span><span class="hljs-punctuation">)</span>   <span class="hljs-comment"># install the whole suite at once</span>

library<span class="hljs-punctuation">(</span>tidyverse<span class="hljs-punctuation">)</span>
<span class="hljs-comment"># After loading, automatically includes: ggplot2, dplyr, tidyr, readr, purrr, stringr, forcats, tibble</span>
</code></pre>
<h2>2. The Pipe Operator %>%</h2>
<p>The pipe passes the result on the left as the first argument of the function on the right, turning "nested calls" into a "pipeline":</p>
<pre><code class="hljs language-r"><span class="hljs-comment"># Nested style (hard to read)</span>
<span class="hljs-built_in">round</span><span class="hljs-punctuation">(</span>mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">sqrt</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Pipe style (read from left to right)</span>
<span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> <span class="hljs-built_in">sqrt</span><span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> mean<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> <span class="hljs-built_in">round</span><span class="hljs-punctuation">(</span><span class="hljs-number">2</span><span class="hljs-punctuation">)</span>
</code></pre>
<ul>
<li>Both <code>%>%</code> (magrittr) and <code>|></code> (native, R 4.1+) are available; the tidyverse ecosystem conventionally uses <code>%>%</code></li>
<li>Shortcut: In RStudio, press <code>Ctrl+Shift+M</code> to insert <code>%>%</code></li>
</ul>
<h2>3. Data Preparation</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># Create example data: gene expression experiment</span>
expr_data <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>
  gene      <span class="hljs-operator">=</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"MYC"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> each <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  condition <span class="hljs-operator">=</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"treatment"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> length.out <span class="hljs-operator">=</span> <span class="hljs-number">12</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  replicate <span class="hljs-operator">=</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  <span class="hljs-built_in">expression</span> <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">10.2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">11.1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">9.8</span><span class="hljs-punctuation">,</span>    <span class="hljs-comment"># TP53 control</span>
                 <span class="hljs-number">12.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">13.0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">12.1</span><span class="hljs-punctuation">,</span>   <span class="hljs-comment"># TP53 treatment</span>
                 <span class="hljs-number">8.1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7.9</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.4</span><span class="hljs-punctuation">,</span>      <span class="hljs-comment"># BRCA1</span>
                 <span class="hljs-number">45.0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">44.2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">45.8</span><span class="hljs-punctuation">,</span>
                 <span class="hljs-number">33.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">32.9</span><span class="hljs-punctuation">,</span> <span class="hljs-number">33.1</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>
print<span class="hljs-punctuation">(</span>expr_data<span class="hljs-punctuation">)</span>
</code></pre>
<p><code>tibble</code> is tidyverse's enhanced data frame: prints more friendlily, and does not automatically convert strings to factors.</p>
<h2>4. Core Verbs</h2>
<h3>4.1 filter: Select Rows</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Samples with expression greater than 20</span>
expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">20</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Multiple conditions: &#x26; for AND, | for OR</span>
expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span>condition <span class="hljs-operator">==</span> <span class="hljs-string">"treatment"</span> <span class="hljs-operator">&#x26;</span> <span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">12</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># %in% for matching a set</span>
expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span>gene <span class="hljs-operator">%in%</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"MYC"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.2 select: Choose Columns</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  select<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> <span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Exclude columns</span>
expr_data <span class="hljs-operator">%>%</span>
  select<span class="hljs-punctuation">(</span><span class="hljs-operator">-</span>replicate<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Helper selectors</span>
expr_data <span class="hljs-operator">%>%</span>
  select<span class="hljs-punctuation">(</span>starts_with<span class="hljs-punctuation">(</span><span class="hljs-string">"expr"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>      <span class="hljs-comment"># column names starting with "expr"</span>
</code></pre>
<h3>4.3 mutate: Add/Modify Columns</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>
    log2expr <span class="hljs-operator">=</span> log2<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    group <span class="hljs-operator">=</span> paste<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> condition<span class="hljs-punctuation">,</span> sep <span class="hljs-operator">=</span> <span class="hljs-string">"_"</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Conditional column generation: case_when</span>
expr_data <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>level <span class="hljs-operator">=</span> case_when<span class="hljs-punctuation">(</span>
    <span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">30</span> <span class="hljs-operator">~</span> <span class="hljs-string">"high"</span><span class="hljs-punctuation">,</span>
    <span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">10</span> <span class="hljs-operator">~</span> <span class="hljs-string">"medium"</span><span class="hljs-punctuation">,</span>
    <span class="hljs-literal">TRUE</span>            <span class="hljs-operator">~</span> <span class="hljs-string">"low"</span>
  <span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.4 arrange: Sort Rows</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>              <span class="hljs-comment"># ascending</span>

expr_data <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>desc<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>        <span class="hljs-comment"># descending</span>

expr_data <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> desc<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>  <span class="hljs-comment"># multi-level sorting</span>
</code></pre>
<h3>4.5 summarise: Summary Statistics</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>
    mean_expr <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    sd_expr <span class="hljs-operator">=</span> sd<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    n <span class="hljs-operator">=</span> n<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.6 group_by + summarise: Grouped Summaries (Most Common Combination)</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> condition<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>
    mean_expr <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    sd_expr <span class="hljs-operator">=</span> sd<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    n <span class="hljs-operator">=</span> n<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span>
  <span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.7 A Full Pipe Chain</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Average expression per gene under treatment, top 3 in descending order</span>
top_genes <span class="hljs-operator">&#x3C;-</span> expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span>condition <span class="hljs-operator">==</span> <span class="hljs-string">"treatment"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>mean_expr <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>desc<span class="hljs-punctuation">(</span>mean_expr<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  slice_max<span class="hljs-punctuation">(</span>mean_expr<span class="hljs-punctuation">,</span> n <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span>

print<span class="hljs-punctuation">(</span>top_genes<span class="hljs-punctuation">)</span>
</code></pre>
<h2>5. Joining Multiple Tables</h2>
<pre><code class="hljs language-r">gene_info <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"MYC"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  chromosome <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  function_desc <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"tumor suppressor"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"DNA repair"</span><span class="hljs-punctuation">,</span>
                    <span class="hljs-string">"receptor kinase"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"transcription factor"</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Left join: keeps all rows from the left table</span>
expr_data <span class="hljs-operator">%>%</span>
  left_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Inner join: keeps only rows found in both tables</span>
expr_data <span class="hljs-operator">%>%</span>
  inner_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Right join / Full join</span>
expr_data <span class="hljs-operator">%>%</span> right_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>
expr_data <span class="hljs-operator">%>%</span> full_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>6. tidyr: Data Reshaping</h2>
<h3>6.1 pivot_longer: Wide to Long</h3>
<pre><code class="hljs language-r">wide <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  sample1 <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">10.2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  sample2 <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">12.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7.9</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  sample3 <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">9.8</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.4</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

long <span class="hljs-operator">&#x3C;-</span> wide <span class="hljs-operator">%>%</span>
  pivot_longer<span class="hljs-punctuation">(</span>cols <span class="hljs-operator">=</span> starts_with<span class="hljs-punctuation">(</span><span class="hljs-string">"sample"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
               names_to <span class="hljs-operator">=</span> <span class="hljs-string">"sample"</span><span class="hljs-punctuation">,</span>
               values_to <span class="hljs-operator">=</span> <span class="hljs-string">"expression"</span><span class="hljs-punctuation">)</span>
print<span class="hljs-punctuation">(</span>long<span class="hljs-punctuation">)</span>
</code></pre>
<h3>6.2 pivot_wider: Long to Wide</h3>
<pre><code class="hljs language-r">long <span class="hljs-operator">%>%</span>
  pivot_wider<span class="hljs-punctuation">(</span>names_from <span class="hljs-operator">=</span> sample<span class="hljs-punctuation">,</span> values_from <span class="hljs-operator">=</span> <span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>6.3 Handling Missing Values</h3>
<pre><code class="hljs language-r">df <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>a <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-literal">NA</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> b <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-literal">NA</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">6</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>

df <span class="hljs-operator">%>%</span> drop_na<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>              <span class="hljs-comment"># remove rows with missing values</span>
df <span class="hljs-operator">%>%</span> replace_na<span class="hljs-punctuation">(</span><span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span>a <span class="hljs-operator">=</span> <span class="hljs-number">0</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>  <span class="hljs-comment"># fill specified columns</span>
df <span class="hljs-operator">%>%</span> fill<span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span>                <span class="hljs-comment"># fill downward</span>
</code></pre>
<h2>7. Practical Example: A Complete Differential Analysis Workflow</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># 1. Read data</span>
<span class="hljs-comment"># expr_data &#x3C;- read_csv("expression.csv")</span>

<span class="hljs-comment"># 2. Data cleaning</span>
cleaned <span class="hljs-operator">&#x3C;-</span> expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span><span class="hljs-operator">!</span><span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>              <span class="hljs-comment"># remove missing values</span>
  mutate<span class="hljs-punctuation">(</span>log2expr <span class="hljs-operator">=</span> log2<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>     <span class="hljs-comment"># transformation</span>
  left_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>           <span class="hljs-comment"># annotation</span>

<span class="hljs-comment"># 3. Summary statistics</span>
summary_table <span class="hljs-operator">&#x3C;-</span> cleaned <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> condition<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>
    mean_log2 <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span>log2expr<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    n <span class="hljs-operator">=</span> n<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span>
  <span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 4. Calculate differences between treatment and control</span>
diff_table <span class="hljs-operator">&#x3C;-</span> summary_table <span class="hljs-operator">%>%</span>
  pivot_wider<span class="hljs-punctuation">(</span>names_from <span class="hljs-operator">=</span> condition<span class="hljs-punctuation">,</span> values_from <span class="hljs-operator">=</span> mean_log2<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>log2fc <span class="hljs-operator">=</span> treatment <span class="hljs-operator">-</span> control<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>desc<span class="hljs-punctuation">(</span>log2fc<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>

print<span class="hljs-punctuation">(</span>diff_table<span class="hljs-punctuation">)</span>
</code></pre>
<h2>8. Summary</h2>
<ul>
<li><strong>Pipes</strong>: <code>%>%</code> lets data flow from left to right, one verb per step</li>
<li><strong>The six key dplyr verbs</strong>: <code>filter</code> / <code>select</code> / <code>mutate</code> / <code>arrange</code> / <code>summarise</code> / <code>group_by</code></li>
<li><strong>Join family</strong>: <code>left_join</code> / <code>inner_join</code>, etc., merge tables by keys</li>
<li><strong>tidyr</strong>: <code>pivot_longer</code> / <code>pivot_wider</code> convert between wide and long tables</li>
</ul>
<p>The next article will introduce ggplot2: creating publication-quality graphics using the grammar of layers.</p>`,qb=`<h1>R tidyverse 数据操作：dplyr 管道与数据清洗</h1>
<p>tidyverse 是一套风格统一的 R 包集合，其中 <code>dplyr</code> 负责数据操作、<code>tidyr</code> 负责数据整形。核心思想是<strong>管道串联动词</strong>，让数据处理的每一步都清晰可读。</p>
<h2>1. 安装与加载</h2>
<pre><code class="hljs language-r">install.packages<span class="hljs-punctuation">(</span><span class="hljs-string">"tidyverse"</span><span class="hljs-punctuation">)</span>   <span class="hljs-comment"># 一次安装全家桶</span>

library<span class="hljs-punctuation">(</span>tidyverse<span class="hljs-punctuation">)</span>
<span class="hljs-comment"># 加载后自动附带: ggplot2, dplyr, tidyr, readr, purrr, stringr, forcats, tibble</span>
</code></pre>
<h2>2. 管道符 %>%</h2>
<p>管道把左边的结果作为右边函数的第一个参数，把"嵌套调用"变成"流水线"：</p>
<pre><code class="hljs language-r"><span class="hljs-comment"># 嵌套写法（难读）</span>
<span class="hljs-built_in">round</span><span class="hljs-punctuation">(</span>mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">sqrt</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 管道写法（从左到右读）</span>
<span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> <span class="hljs-built_in">sqrt</span><span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> mean<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> <span class="hljs-built_in">round</span><span class="hljs-punctuation">(</span><span class="hljs-number">2</span><span class="hljs-punctuation">)</span>
</code></pre>
<ul>
<li><code>%>%</code>（magrittr）与 <code>|></code>（R 4.1+ 原生）都可用；tidyverse 生态惯用 <code>%>%</code></li>
<li>快捷键：RStudio 中 <code>Ctrl+Shift+M</code> 插入 <code>%>%</code></li>
</ul>
<h2>3. 数据准备</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># 创建示例数据：基因表达实验</span>
expr_data <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>
  gene      <span class="hljs-operator">=</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"MYC"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> each <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  condition <span class="hljs-operator">=</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"treatment"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> length.out <span class="hljs-operator">=</span> <span class="hljs-number">12</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  replicate <span class="hljs-operator">=</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  <span class="hljs-built_in">expression</span> <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">10.2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">11.1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">9.8</span><span class="hljs-punctuation">,</span>    <span class="hljs-comment"># TP53 control</span>
                 <span class="hljs-number">12.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">13.0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">12.1</span><span class="hljs-punctuation">,</span>   <span class="hljs-comment"># TP53 treatment</span>
                 <span class="hljs-number">8.1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7.9</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.4</span><span class="hljs-punctuation">,</span>      <span class="hljs-comment"># BRCA1</span>
                 <span class="hljs-number">45.0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">44.2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">45.8</span><span class="hljs-punctuation">,</span>
                 <span class="hljs-number">33.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">32.9</span><span class="hljs-punctuation">,</span> <span class="hljs-number">33.1</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>
print<span class="hljs-punctuation">(</span>expr_data<span class="hljs-punctuation">)</span>
</code></pre>
<p><code>tibble</code> 是 tidyverse 的增强数据框：打印更友好、不自动转字符串为因子。</p>
<h2>4. 核心动词</h2>
<h3>4.1 filter：筛选行</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 表达量大于 20 的样本</span>
expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">20</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 多条件：&#x26; 且，| 或</span>
expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span>condition <span class="hljs-operator">==</span> <span class="hljs-string">"treatment"</span> <span class="hljs-operator">&#x26;</span> <span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">12</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># %in% 匹配集合</span>
expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span>gene <span class="hljs-operator">%in%</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"MYC"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.2 select：选择列</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  select<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> <span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 排除列</span>
expr_data <span class="hljs-operator">%>%</span>
  select<span class="hljs-punctuation">(</span><span class="hljs-operator">-</span>replicate<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 辅助选择器</span>
expr_data <span class="hljs-operator">%>%</span>
  select<span class="hljs-punctuation">(</span>starts_with<span class="hljs-punctuation">(</span><span class="hljs-string">"expr"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>      <span class="hljs-comment"># 列名以 expr 开头</span>
</code></pre>
<h3>4.3 mutate：新增/修改列</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>
    log2expr <span class="hljs-operator">=</span> log2<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    group <span class="hljs-operator">=</span> paste<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> condition<span class="hljs-punctuation">,</span> sep <span class="hljs-operator">=</span> <span class="hljs-string">"_"</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 条件生成列：case_when</span>
expr_data <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>level <span class="hljs-operator">=</span> case_when<span class="hljs-punctuation">(</span>
    <span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">30</span> <span class="hljs-operator">~</span> <span class="hljs-string">"high"</span><span class="hljs-punctuation">,</span>
    <span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">10</span> <span class="hljs-operator">~</span> <span class="hljs-string">"medium"</span><span class="hljs-punctuation">,</span>
    <span class="hljs-literal">TRUE</span>            <span class="hljs-operator">~</span> <span class="hljs-string">"low"</span>
  <span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.4 arrange：排序</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>              <span class="hljs-comment"># 升序</span>

expr_data <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>desc<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>        <span class="hljs-comment"># 降序</span>

expr_data <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> desc<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>  <span class="hljs-comment"># 多级排序</span>
</code></pre>
<h3>4.5 summarise：汇总</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>
    mean_expr <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    sd_expr <span class="hljs-operator">=</span> sd<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    n <span class="hljs-operator">=</span> n<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.6 group_by + summarise：分组汇总（最常用组合）</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> condition<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>
    mean_expr <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    sd_expr <span class="hljs-operator">=</span> sd<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    n <span class="hljs-operator">=</span> n<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span>
  <span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.7 完整管道链</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 每个基因在 treatment 下的平均表达，按降序取前 3</span>
top_genes <span class="hljs-operator">&#x3C;-</span> expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span>condition <span class="hljs-operator">==</span> <span class="hljs-string">"treatment"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>mean_expr <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>desc<span class="hljs-punctuation">(</span>mean_expr<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  slice_max<span class="hljs-punctuation">(</span>mean_expr<span class="hljs-punctuation">,</span> n <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span>

print<span class="hljs-punctuation">(</span>top_genes<span class="hljs-punctuation">)</span>
</code></pre>
<h2>5. 连接多个表（join）</h2>
<pre><code class="hljs language-r">gene_info <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"MYC"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  chromosome <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  function_desc <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"tumor suppressor"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"DNA repair"</span><span class="hljs-punctuation">,</span>
                    <span class="hljs-string">"receptor kinase"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"transcription factor"</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 左连接：保留左表所有行</span>
expr_data <span class="hljs-operator">%>%</span>
  left_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 内连接：只保留两边都有的</span>
expr_data <span class="hljs-operator">%>%</span>
  inner_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 右连接 / 全连接</span>
expr_data <span class="hljs-operator">%>%</span> right_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>
expr_data <span class="hljs-operator">%>%</span> full_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>6. tidyr：数据整形</h2>
<h3>6.1 pivot_longer：宽表转长表</h3>
<pre><code class="hljs language-r">wide <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  sample1 <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">10.2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  sample2 <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">12.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7.9</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  sample3 <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">9.8</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.4</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

long <span class="hljs-operator">&#x3C;-</span> wide <span class="hljs-operator">%>%</span>
  pivot_longer<span class="hljs-punctuation">(</span>cols <span class="hljs-operator">=</span> starts_with<span class="hljs-punctuation">(</span><span class="hljs-string">"sample"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
               names_to <span class="hljs-operator">=</span> <span class="hljs-string">"sample"</span><span class="hljs-punctuation">,</span>
               values_to <span class="hljs-operator">=</span> <span class="hljs-string">"expression"</span><span class="hljs-punctuation">)</span>
print<span class="hljs-punctuation">(</span>long<span class="hljs-punctuation">)</span>
</code></pre>
<h3>6.2 pivot_wider：长表转宽表</h3>
<pre><code class="hljs language-r">long <span class="hljs-operator">%>%</span>
  pivot_wider<span class="hljs-punctuation">(</span>names_from <span class="hljs-operator">=</span> sample<span class="hljs-punctuation">,</span> values_from <span class="hljs-operator">=</span> <span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>6.3 缺失值处理</h3>
<pre><code class="hljs language-r">df <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>a <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-literal">NA</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> b <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-literal">NA</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">6</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>

df <span class="hljs-operator">%>%</span> drop_na<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>              <span class="hljs-comment"># 删除含缺失的行</span>
df <span class="hljs-operator">%>%</span> replace_na<span class="hljs-punctuation">(</span><span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span>a <span class="hljs-operator">=</span> <span class="hljs-number">0</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>  <span class="hljs-comment"># 指定列填充</span>
df <span class="hljs-operator">%>%</span> fill<span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span>                <span class="hljs-comment"># 向下填充</span>
</code></pre>
<h2>7. 实战：完整差异分析工作流</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># 1. 读取数据</span>
<span class="hljs-comment"># expr_data &#x3C;- read_csv("expression.csv")</span>

<span class="hljs-comment"># 2. 数据清洗</span>
cleaned <span class="hljs-operator">&#x3C;-</span> expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span><span class="hljs-operator">!</span><span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>              <span class="hljs-comment"># 去缺失</span>
  mutate<span class="hljs-punctuation">(</span>log2expr <span class="hljs-operator">=</span> log2<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>     <span class="hljs-comment"># 转换</span>
  left_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>           <span class="hljs-comment"># 注释</span>

<span class="hljs-comment"># 3. 汇总统计</span>
summary_table <span class="hljs-operator">&#x3C;-</span> cleaned <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> condition<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>
    mean_log2 <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span>log2expr<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    n <span class="hljs-operator">=</span> n<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span>
  <span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 4. 计算处理/对照的差异</span>
diff_table <span class="hljs-operator">&#x3C;-</span> summary_table <span class="hljs-operator">%>%</span>
  pivot_wider<span class="hljs-punctuation">(</span>names_from <span class="hljs-operator">=</span> condition<span class="hljs-punctuation">,</span> values_from <span class="hljs-operator">=</span> mean_log2<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>log2fc <span class="hljs-operator">=</span> treatment <span class="hljs-operator">-</span> control<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>desc<span class="hljs-punctuation">(</span>log2fc<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>

print<span class="hljs-punctuation">(</span>diff_table<span class="hljs-punctuation">)</span>
</code></pre>
<h2>8. 小结</h2>
<ul>
<li><strong>管道</strong>：<code>%>%</code> 让数据流从左到右，每一步一个动词</li>
<li><strong>dplyr 六大动词</strong>：<code>filter</code> / <code>select</code> / <code>mutate</code> / <code>arrange</code> / <code>summarise</code> / <code>group_by</code></li>
<li><strong>join 家族</strong>：<code>left_join</code> / <code>inner_join</code> 等按键合并表</li>
<li><strong>tidyr</strong>：<code>pivot_longer</code> / <code>pivot_wider</code> 在宽长表之间转换</li>
</ul>
<p>下一篇将介绍 ggplot2：用图层语法绘制出版级图表。</p>`,zb=`<h1>Review of Cryo-EM Technology</h1>
<p>In 2015, Nature Methods named cryo-electron microscopy (cryo-EM) its "Method of the Year." A decade on, this technique has evolved from the subject of ridicule as "blobology" into a mainstay for determining the structures of biological macromolecules. This article reviews the key milestones of this technological revolution and looks ahead to structural biology in the AI era.</p>
<h2>1. The Eve of Revolution: Why Cryo-EM Was Once Called "Blobology"</h2>
<h3>1.1 The Three Pillars of Traditional Structural Biology</h3>
<ul>
<li><strong>X-ray crystallography</strong>: Requires high-quality crystals; many membrane proteins and large complexes are difficult to crystallize</li>
<li><strong>Nuclear magnetic resonance (NMR)</strong>: Limited by molecular weight (typically &#x3C; 50 kDa)</li>
<li><strong>Cryo-EM</strong>: No upper molecular weight limit in theory, no crystals required, but resolution remained stuck at 10–30 Å for a long time</li>
</ul>
<h3>1.2 Root Causes of the Resolution Bottleneck</h3>
<ol>
<li><strong>Electron damage</strong>: The electron beam damages the sample, limiting the total dose to ~20 e⁻/Å², resulting in extremely low signal-to-noise ratio</li>
<li><strong>Outdated detectors</strong>: CCD cameras are noisy at low doses and lack single-electron counting capability</li>
<li><strong>Specimen drift</strong>: Electron-beam-induced specimen movement blurs the images</li>
<li><strong>Missing phase information</strong>: Under the weak-phase-object approximation, defocused imaging loses part of the phase information</li>
</ol>
<h2>2. The Three Cornerstones of the Technological Revolution</h2>
<h3>2.1 Direct Electron Detectors (DET) — 2013–2015</h3>
<table>
<thead>
<tr>
<th>Detector</th>
<th>Year Released</th>
<th>Key Features</th>
</tr>
</thead>
<tbody>
<tr>
<td>Gatan K2 Summit</td>
<td>2013</td>
<td>First commercial CMOS direct detector, single-electron counting</td>
</tr>
<tr>
<td>FEI Falcon II</td>
<td>2014</td>
<td>High frame rate, lossless readout</td>
</tr>
<tr>
<td>Gatan K3</td>
<td>2017</td>
<td>Larger field of view, super-resolution</td>
</tr>
<tr>
<td>Falcon 4 / Falcon 4i</td>
<td>2020</td>
<td>Electron counting + integrated energy filtering</td>
</tr>
</tbody>
</table>
<p><strong>Why DET is revolutionary</strong>:</p>
<ul>
<li>Single-electron counting mode eliminates readout noise → greatly improved signal-to-noise ratio</li>
<li>Ultra-fast frame rates (hundreds of frames per second) → <strong>movie mode</strong>, enabling per-frame correction of specimen drift</li>
<li>Combined with energy filters (GIF/Selectris), removes inelastic scattering background</li>
</ul>
<h3>2.2 The Stability Revolution — 2016–2018</h3>
<p>The "Resolution Revolution" was not just about detectors; it was equally a triumph of <strong>overall stability</strong>:</p>
<ul>
<li>Constant temperature and humidity in the microscope room (±0.1 °C), soundproofing and vibration isolation</li>
<li>Improved specimen holders and grid-clamping mechanisms (e.g., autoloader)</li>
<li>Widespread adoption of dose-symmetric acquisition schemes and low-dose modes</li>
<li>Software: MotionCor2's <strong>per-frame drift correction</strong> preserved high-frequency signals</li>
</ul>
<p>Representative achievements (2017–2020):</p>
<ul>
<li>Human γ-secretase complex (3.4 Å, 2015, Scheres group)</li>
<li>β-galactosidase at 2.2 Å (2016)</li>
<li><strong>Human ribosome at ~2.5 Å</strong> (multiple groups)</li>
<li>Multiple membrane protein channels and receptor complexes broke the 3 Å barrier</li>
</ul>
<h3>2.3 Algorithms and Computing Power — 2017–2021</h3>
<ul>
<li><strong>Bayesian methods</strong> (RELION) made 3D classification/refinement robust</li>
<li><strong>cryoSPARC's ab-initio and NU-refinement</strong>: heterogeneous species separation, flexible-region refinement</li>
<li><strong>GPU acceleration</strong> made iterative reconstruction of million-particle datasets routine</li>
<li>Deep learning enters: <strong>Topaz / crYOLO</strong> for particle picking, <strong>DeepEM</strong> for denoising</li>
</ul>
<h2>3. The 2020s: From 3 Å to Atomic Resolution</h2>
<h3>3.1 Breaking the 2 Å Barrier</h3>
<ul>
<li>2020–2021: Multiple samples reached 1.5–2.0 Å (e.g., apoferritin at 1.22 Å, 2020)</li>
<li>Sub-2 Å density maps reveal: hydrogen atoms, water molecule networks, side-chain microenvironments, ligand interaction details</li>
<li>This paves the way for <strong>density-based de novo modeling</strong> and drug design</li>
</ul>
<h3>3.2 Time-Resolved Cryo-EM</h3>
<ul>
<li>Mix-spray techniques capture enzyme catalytic intermediates</li>
<li>Microfluidics and rapid freezing (&#x3C; 10 ms timescale)</li>
<li>Goal: observe conformational change dynamics, rather than a single static structure</li>
</ul>
<h3>3.3 In Situ Structural Biology (In situ / Cellular Cryo-EM)</h3>
<ul>
<li>Cryo-focused ion beam (cryo-FIB) thinning + cryo-electron tomography (cryo-ET)</li>
<li>Resolving the assembly states and spatial distribution of protein complexes <strong>within cells</strong></li>
<li>Representative: in situ structures of ribosomes, proteasomes, and vesicle trafficking complexes in cells</li>
</ul>
<h2>4. Synergy and Challenges in the AI Era</h2>
<h3>4.1 The AlphaFold Shockwave (2021–)</h3>
<p>AlphaFold2 / AlphaFold3 can predict protein structures with high accuracy, forming a <strong>complementary rather than replacement</strong> relationship with cryo-EM:</p>
<table>
<thead>
<tr>
<th>Dimension</th>
<th>AlphaFold</th>
<th>Cryo-EM</th>
</tr>
</thead>
<tbody>
<tr>
<td>Input</td>
<td>Sequence</td>
<td>Experimental sample</td>
</tr>
<tr>
<td>Output</td>
<td>Predicted model</td>
<td>Experimental density/model</td>
</tr>
<tr>
<td>Conformation</td>
<td>Single (or few)</td>
<td>True conformational distribution</td>
</tr>
<tr>
<td>Modifications/ligands</td>
<td>Limited predictive capability</td>
<td>Direct observation</td>
</tr>
<tr>
<td>Complex assembly</td>
<td>Predicted</td>
<td>True in situ state</td>
</tr>
</tbody>
</table>
<p><strong>Actual workflow</strong>: AF2 predicted model → used as initial template for molecular replacement/model building in cryo-EM → refined against experimental density → yielding a model containing true interactions and modifications.</p>
<h3>4.2 Further Penetration of Deep Learning</h3>
<ul>
<li>Automatic density map model building (ModelAngelo, DeepTracer)</li>
<li>End-to-end models for particle picking and denoising</li>
<li>Diffusion-model-based density map deblurring</li>
<li><strong>Joint refinement of AlphaFold and density maps</strong> (e.g., ISOLDE + AF constraints)</li>
</ul>
<h3>4.3 Challenges Remain</h3>
<ul>
<li><strong>Small proteins (&#x3C; 50 kDa)</strong>: low signal-to-noise ratio, difficult to improve resolution (require better phase plates/detectors)</li>
<li><strong>Flexible large complexes</strong>: conformational continua are difficult to describe with discrete classification (3DVA, 3D FLEX, etc. are addressing this)</li>
<li><strong>Membrane proteins</strong>: detergent/nanodisc environments differ from native membranes</li>
<li><strong>Specimen preparation</strong>: the air-water interface problem remains one of the major bottlenecks</li>
</ul>
<h2>5. Future Outlook</h2>
<ol>
<li><strong>Resolution becoming routine</strong>: sub-2 Å transitions from "world record" to a routine goal</li>
<li><strong>Dynamic structural biology</strong>: time-resolved + continuous conformational analysis, moving from "snapshots" to "movies"</li>
<li><strong>In-cell structural biology</strong>: cryo-ET combined with FIB to understand molecular machines in the cellular environment</li>
<li><strong>Full AI integration</strong>: automated, intelligent pipelines from picking to modeling</li>
<li><strong>Deep integration with drug design</strong>: high-resolution structures + virtual screening + generative design closed loop</li>
</ol>
<h2>6. Summary</h2>
<ul>
<li>Direct electron detectors, the stability revolution, and algorithmic advances are the three cornerstones of the "Resolution Revolution"</li>
<li>Cryo-EM has transformed from a "last resort" into the <strong>preferred method</strong> for membrane proteins and large complexes</li>
<li>AlphaFold and cryo-EM are complementary: prediction guides experiment, experiment corrects prediction</li>
<li>Future directions: dynamics, in situ, atomic resolution, AI integration</li>
</ul>
<p>The past decade of cryo-EM is a model of synergy among engineering, physics, and computational biology. In the next decade, it will join forces with AI to redefine how we understand life.</p>`,Vb=`<h1>冷冻电镜技术综述</h1>
<p>2015 年，Nature Methods 将冷冻电镜（cryo-EM）评为"年度方法"。十年过去，这项技术从"blobology"（模糊团块学）的调侃对象，成长为解析生物大分子结构的主力军。本文回顾这场技术革命的关键节点，并展望 AI 时代的结构生物学。</p>
<h2>1. 革命前夜：为什么冷冻电镜曾被称为"blobology"</h2>
<h3>1.1 传统结构生物学的三驾马车</h3>
<ul>
<li><strong>X 射线晶体学</strong>：需要高质量晶体，许多膜蛋白、大复合物难以结晶</li>
<li><strong>核磁共振（NMR）</strong>：受分子量限制（通常 &#x3C; 50 kDa）</li>
<li><strong>冷冻电镜</strong>：理论上无分子量上限、无需晶体，但分辨率长期卡在 10–30 Å</li>
</ul>
<h3>1.2 分辨率瓶颈的根源</h3>
<ol>
<li><strong>电子损伤</strong>：电子束会破坏样品，总剂量被限制在 ~20 e⁻/Å²，信噪比极低</li>
<li><strong>探测器落后</strong>：CCD 相机在低剂量下噪声大，且无单电子计数能力</li>
<li><strong>样品漂移</strong>：电子束导致的样品移动模糊图像</li>
<li><strong>相位信息缺失</strong>：弱相位物体近似下，欠焦成像丢失部分相位信息</li>
</ol>
<h2>2. 技术革命的三块基石</h2>
<h3>2.1 直接电子探测器（DET）——2013–2015</h3>
<table>
<thead>
<tr>
<th>探测器</th>
<th>推出年份</th>
<th>关键特性</th>
</tr>
</thead>
<tbody>
<tr>
<td>Gatan K2 Summit</td>
<td>2013</td>
<td>首款商用 CMOS 直接探测，单电子计数</td>
</tr>
<tr>
<td>FEI Falcon II</td>
<td>2014</td>
<td>高帧率、无损读出</td>
</tr>
<tr>
<td>Gatan K3</td>
<td>2017</td>
<td>更大视野、超分辨率</td>
</tr>
<tr>
<td>Falcon 4 / Falcon 4i</td>
<td>2020</td>
<td>电子计数 + 能量过滤集成</td>
</tr>
</tbody>
</table>
<p><strong>为什么 DET 是革命性的</strong>：</p>
<ul>
<li>单电子计数模式消除读出噪声 → 信噪比大幅提升</li>
<li>超快帧率（每秒数百帧）→ <strong>movie 模式</strong>，可逐帧校正样品漂移</li>
<li>与能量过滤器（GIF/Selectris）结合，去除非弹性散射本底</li>
</ul>
<h3>2.2 稳定性革命——2016–2018</h3>
<p>"Resolution Revolution" 不只是探测器，更是<strong>整体稳定性</strong>的胜利：</p>
<ul>
<li>电镜室恒温恒湿（±0.1 °C）、隔音防震</li>
<li>样品杆与载网夹持机构改进（如 autoloader）</li>
<li>剂量对称采集方案、低剂量模式的普及</li>
<li>软件层面：MotionCor2 的<strong>逐帧漂移校正</strong>让高频信号得以保留</li>
</ul>
<p>代表性成果（2017–2020）：</p>
<ul>
<li>人源 γ-分泌酶复合物（3.4 Å，2015，Scheres 组）</li>
<li>β-半乳糖苷酶 2.2 Å（2016）</li>
<li><strong>人源核糖体 2.5 Å 级</strong>（多组）</li>
<li>多种膜蛋白通道、受体复合物突破 3 Å</li>
</ul>
<h3>2.3 算法与算力——2017–2021</h3>
<ul>
<li><strong>贝叶斯方法</strong>（RELION）让 3D 分类/精修稳健化</li>
<li><strong>cryoSPARC 的 Ab-initio 与 NU-refinement</strong>：异构体分离、柔性区域精修</li>
<li><strong>GPU 加速</strong>使百万级颗粒的迭代重构成为常规操作</li>
<li>深度学习进入：<strong>Topaz / crYOLO</strong> 颗粒挑选、<strong>DeepEM</strong> 去噪</li>
</ul>
<h2>3. 2020 年代：从 3 Å 到原子级</h2>
<h3>3.1 打破 2 Å 门槛</h3>
<ul>
<li>2020–2021：多种样品达到 1.5–2.0 Å（如去铁蛋白 1.22 Å，2020）</li>
<li>亚 2 Å 密度下可观察到：氢原子、水分子网络、侧链微环境、配体互作细节</li>
<li>这为<strong>基于密度的从头建模</strong>与药物设计铺平道路</li>
</ul>
<h3>3.2 时间分辨冷冻电镜（Time-resolved cryo-EM）</h3>
<ul>
<li>混合-喷雾（mix-spray）技术捕捉酶催化中间态</li>
<li>微流控与快速冷冻（&#x3C; 10 ms 时间尺度）</li>
<li>目标：观察构象变化动态过程，而非单一静态结构</li>
</ul>
<h3>3.3 原位结构生物学（In situ / Cellular cryo-EM）</h3>
<ul>
<li>冷冻聚焦离子束（cryo-FIB）减薄 + 冷冻电子断层扫描（cryo-ET）</li>
<li>在<strong>细胞原位</strong>解析蛋白质复合物的组装状态与空间分布</li>
<li>代表：细胞内核糖体、蛋白酶体、囊泡运输复合物的原位结构</li>
</ul>
<h2>4. AI 时代的协同与挑战</h2>
<h3>4.1 AlphaFold 冲击波（2021–）</h3>
<p>AlphaFold2 / AlphaFold3 可以高精度预测蛋白结构，与冷冻电镜形成<strong>互补而非替代</strong>：</p>
<table>
<thead>
<tr>
<th>维度</th>
<th>AlphaFold</th>
<th>冷冻电镜</th>
</tr>
</thead>
<tbody>
<tr>
<td>输入</td>
<td>序列</td>
<td>实验样品</td>
</tr>
<tr>
<td>输出</td>
<td>预测模型</td>
<td>实验密度/模型</td>
</tr>
<tr>
<td>构象</td>
<td>单一（或少量）</td>
<td>真实构象分布</td>
</tr>
<tr>
<td>修饰/配体</td>
<td>预测能力有限</td>
<td>直接观察</td>
</tr>
<tr>
<td>复合物组装</td>
<td>预测</td>
<td>原位真实状态</td>
</tr>
</tbody>
</table>
<p><strong>实际工作流</strong>：AF2 预测模型 → 作为冷冻电镜分子替换（molecular replacement）/建模的初始模板 → 实验密度修正 → 获得包含真实互作与修饰的模型。</p>
<h3>4.2 深度学习的进一步渗透</h3>
<ul>
<li>密度图自动建模（ModelAngelo、DeepTracer）</li>
<li>颗粒挑选与去噪的端到端模型</li>
<li>基于扩散模型的密度图去模糊</li>
<li><strong>AlphaFold 与密度图的联合精修</strong>（如 ISOLDE + AF 约束）</li>
</ul>
<h3>4.3 挑战依然存在</h3>
<ul>
<li><strong>小蛋白（&#x3C; 50 kDa）</strong>：信噪比低，分辨率提升困难（需要更好的相板/探测器）</li>
<li><strong>柔性大复合物</strong>：构象连续体难以用离散分类描述（3DVA、3D FLEX 等正在解决）</li>
<li><strong>膜蛋白</strong>：去垢剂/纳米盘环境与天然膜差异</li>
<li><strong>样品制备</strong>：air-water interface 问题仍是主要瓶颈之一</li>
</ul>
<h2>5. 未来图景</h2>
<ol>
<li><strong>分辨率常态化</strong>：亚 2 Å 从"世界纪录"变为常规目标</li>
<li><strong>动态结构生物学</strong>：时间分辨 + 连续构象分析，从"快照"走向"电影"</li>
<li><strong>细胞内结构生物学</strong>：cryo-ET 与 FIB 结合，在细胞环境中理解分子机器</li>
<li><strong>AI 全流程集成</strong>：从挑选到建模的自动化、智能化管线</li>
<li><strong>与药物设计深度融合</strong>：高分辨结构 + 虚拟筛选 + 生成式设计闭环</li>
</ol>
<h2>6. 小结</h2>
<ul>
<li>直接电子探测器、稳定性革命、算法进步是"分辨率革命"的三块基石</li>
<li>冷冻电镜已从"最后手段"变为膜蛋白、大复合物的<strong>首选方法</strong></li>
<li>AlphaFold 与冷冻电镜形成互补：预测引导实验，实验修正预测</li>
<li>未来方向：动态、原位、原子级、AI 融合</li>
</ul>
<p>冷冻电镜的十年，是工程、物理与计算生物学协同的典范。下一个十年，它将与 AI 一起重新定义我们理解生命的方式。</p>`,Ub=`<h1>Cryo-EM Single Particle Analysis: Full Data Processing Workflow</h1>
<p>Cryo-EM single particle analysis (SPA) is a core technique for determining near-atomic resolution structures of biological macromolecules. This article outlines the complete pipeline from micrographs to atomic models, introducing the purpose, common software, and key parameters for each step.</p>
<h2>1. Workflow Overview</h2>
<pre><code>Sample preparation (cryo-grid preparation)
   ↓
Data acquisition (Titan Krios / Glacios, etc.)
   ↓
① Preprocessing: motion correction + CTF estimation
   ↓
② Particle picking
   ↓
③ 2D classification (removing bad particles)
   ↓
④ Initial model (Ab initio)
   ↓
⑤ 3D classification (heterogeneity analysis)
   ↓
⑥ 3D refinement
   ↓
⑦ Resolution assessment (FSC curve)
   ↓
⑧ Atomic model building and refinement
</code></pre>
<h2>2. Sample Preparation and Data Acquisition</h2>
<h3>2.1 Cryo-grid Preparation</h3>
<ul>
<li>Apply the protein sample (typically 0.5–5 mg/mL) to a holey carbon grid</li>
<li>Use a Vitrobot / Leica GP or similar device to rapidly plunge into liquid ethane, forming <strong>vitreous ice</strong></li>
<li>Ideal ice thickness: slightly thinner than the particle diameter; ice that is too thick reduces signal-to-noise, while ice that is too thin may deform particles</li>
</ul>
<h3>2.2 Data Acquisition</h3>
<ul>
<li>300 kV electron microscope (Titan Krios) + direct electron detector (Gatan K3 / Falcon 4)</li>
<li><strong>Super-resolution mode</strong> + dose fractionation (total dose about 40–60 e⁻/Å²)</li>
<li>Each movie typically contains 30–50 frames for motion correction</li>
<li>Target: 10,000–20,000 micrographs (high-resolution structures usually require millions of particles)</li>
</ul>
<h2>3. Preprocessing</h2>
<h3>3.1 Motion Correction</h3>
<p>Sample drift caused by the electron beam blurs images and must be corrected frame by frame:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Motion correction in RELION 4</span>
relion_run_motioncorr \\
  --i Micrographs/ \\
  --o MotionCorr/ \\
  --j 8 \\
  --pipeline_control MotionCorr/job001/
</code></pre>
<p>Common software: <strong>MotionCor2</strong> (classic), <strong>built-in RELION</strong>, <strong>cryoSPARC Patch Motion Correction</strong>.</p>
<h3>3.2 CTF Estimation</h3>
<p>The CTF (Contrast Transfer Function) describes the phase flipping and defocus effects of electron microscope imaging:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># CTFFIND4 (classic command line)</span>
ctffind --input_micrograph MotionCorr/mic1.mrc \\
        --output_diag mic1_diag.mrc \\
        --defU 20000 --defV 20000

<span class="hljs-comment"># Gctf (faster GPU version)</span>
Gctf --input <span class="hljs-string">'MotionCorr/*.mrc'</span> --apix 1.0
</code></pre>
<p>Common software: <strong>CTFFIND4</strong>, <strong>Gctf</strong>, <strong>cryoSPARC Patch CTF</strong>.</p>
<p><strong>Key criteria</strong>: The resolution cutoff is usually taken as 3–4 Å; discard micrographs with excessive astigmatism or severe drift.</p>
<h2>4. Particle Picking</h2>
<h3>4.1 Traditional Methods</h3>
<ul>
<li><strong>Template-based</strong>: use 2D class averages as templates and perform correlation searches</li>
<li>Software: auto-picking in RELION, Template Picker in cryoSPARC</li>
</ul>
<h3>4.2 Deep Learning Approaches (currently mainstream)</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Topaz: CNN-based particle picking</span>
topaz train --model-path model.pkl --epochs 20 \\
  --num-workers 4 --train-data particles.toml

topaz extract --model-path model.pkl \\
  --micrographs micrographs.toml \\
  --output-prefix picked --radius 100
</code></pre>
<p>Mainstream tools: <strong>Topaz</strong>, <strong>crYOLO</strong> (very fast, supports real-time picking), <strong>cryoSPARC Blob Picker</strong>.</p>
<blockquote>
<p>In modern workflows, a blob/template picker is often used first to quickly select a subset for 2D classification, and clean class averages are then used to train Topaz/crYOLO models, significantly improving picking quality for complex samples (membrane proteins, small proteins).</p>
</blockquote>
<h2>5. 2D Classification: Particle Screening</h2>
<p>2D classification clusters particles by projection orientation while <strong>removing ice contamination, mis-picked particles, and denatured particles</strong>:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RELION 2D classification</span>
relion_refine --o Class2D/ \\
  --particle_images particles.star \\
  --ctf --K 100 --iter 25 --tau2_fudge 4 \\
  --flatten_solvent --zero_mask --pool 30 \\
  --per_particle_ctf --j 8
</code></pre>
<p><strong>Quality criteria</strong>:</p>
<ul>
<li>Class averages should show clear secondary structure (α-helices, β-sheets)</li>
<li>Retain 50–80% of the total particles as "good"</li>
<li>High-quality class averages are a prerequisite for subsequent 3D reconstruction</li>
</ul>
<h2>6. 3D Reconstruction: Initial Model and Classification</h2>
<h3>6.1 Initial Model (Ab initio)</h3>
<p>No reference model is required; reconstruction starts from random initial structures:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># cryoSPARC</span>
cryosparc_abinitio ... --num_classes 3
</code></pre>
<ul>
<li>Commonly used: <strong>cryoSPARC Ab-initio</strong> (fast and robust), RELION 3D initial model</li>
<li>Generate 2–4 classes and select the class with the clearest structural features as a reference</li>
</ul>
<h3>6.2 3D Classification: Handling Conformational Heterogeneity</h3>
<p>Biological samples often exist in multiple conformations (open/closed, ligand-bound/unbound):</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RELION 3D classification</span>
relion_refine --o Class3D/ \\
  --ref initial.mrc --K 4 --iter 40 \\
  --tau2_fudge 8 --particle_diameter 200
</code></pre>
<p><strong>Strategy</strong>: First use coarse classification (K=4–6) to inspect the conformational distribution, then sub-classify the conformation of interest (K=2–3). Ultimately, each 3D class corresponds to one conformation.</p>
<h2>7. 3D Refinement</h2>
<h3>7.1 Standard Refinement</h3>
<pre><code class="hljs language-bash">relion_refine --o Refine3D/ \\
  --ref class3d.mrc --particle_diameter 200 \\
  --auto_mask --solvent_mask \\
  --ctf --per_particle_ctf --j 8
</code></pre>
<p>Key parameters:</p>
<ul>
<li>Start with <strong>C1 symmetry</strong>; for symmetric molecules, specify C/D/O symmetry (e.g., C3, D7)</li>
<li>A <strong>solvent mask</strong> must be used; otherwise the FSC is artificially high</li>
<li>Use local refinement to handle flexible regions</li>
</ul>
<h3>7.2 Bayesian Polishing</h3>
<p>Performs trajectory fitting and weighting at the movie-frame level for particles, improving high-frequency information:</p>
<pre><code class="hljs language-bash">relion_polish --i Refine3D/run_data.star \\
  --o Polish/ \\
  --model_type 2 --nr_iter 5
</code></pre>
<h3>7.3 CTF Refinement</h3>
<p>Refines per-particle defocus, astigmatism, and beam tilt:</p>
<pre><code class="hljs language-bash">relion_ctf_refine --i Polish/particles.star \\
  --o CtfRefine/ \\
  --fit_defocus --fit_magnification --fit_tilt
</code></pre>
<h2>8. Resolution Assessment: FSC Curve</h2>
<p><strong>FSC (Fourier Shell Correlation)</strong>: Particles are randomly split into two half-sets (half-maps), reconstructed separately, and the correlation coefficient between the two maps is computed:</p>
<ul>
<li><strong>FSC = 0.143 criterion</strong>: gold-standard resolution determination (after half-map correction)</li>
<li><strong>FSC = 0.5</strong>: conservative estimate (uncorrected)</li>
</ul>
<pre><code class="hljs language-bash">relion_postprocess --i Refine3D/run_half1_class001.mrc \\
  --mask mask.mrc --autob_masking \\
  --angpix 1.0
</code></pre>
<p><strong>Quality checklist</strong>:</p>
<table>
<thead>
<tr>
<th>Metric</th>
<th>Good standard</th>
</tr>
</thead>
<tbody>
<tr>
<td>Global resolution (FSC 0.143)</td>
<td>≤ 3.5 Å (high resolution)</td>
</tr>
<tr>
<td>Number of particles</td>
<td>≥ 100,000 (hundreds of thousands typically needed for 3 Å)</td>
</tr>
<tr>
<td>Orientation distribution</td>
<td>Uniform (no preferred orientation)</td>
</tr>
<tr>
<td>Density continuity</td>
<td>Side-chain density clearly visible (&#x3C; 3 Å)</td>
</tr>
</tbody>
</table>
<h2>9. Comparison of Common Software Stacks</h2>
<table>
<thead>
<tr>
<th>Step</th>
<th>RELION (command line)</th>
<th>cryoSPARC (GUI)</th>
</tr>
</thead>
<tbody>
<tr>
<td>Motion correction</td>
<td>Built-in / MotionCor2</td>
<td>Patch Motion</td>
</tr>
<tr>
<td>CTF</td>
<td>CTFFIND4 / Gctf</td>
<td>Patch CTF</td>
</tr>
<tr>
<td>Picking</td>
<td>Template / Topaz integration</td>
<td>Blob / Template / Topaz</td>
</tr>
<tr>
<td>2D/3D</td>
<td>Stable and reliable</td>
<td>Fast, interactive</td>
</tr>
<tr>
<td>Non-uniform refinement</td>
<td>—</td>
<td><strong>NU-refinement</strong> (powerful for flexible macromolecules)</td>
</tr>
</tbody>
</table>
<blockquote>
<p><strong>Workflow recommendation</strong>: Use cryoSPARC for preprocessing and initial reconstruction, and RELION for high-precision refinement and polishing; alternatively, use cryoSPARC for the entire workflow (more user-friendly for beginners).</p>
</blockquote>
<h2>10. Summary</h2>
<ul>
<li>Preprocessing (motion correction + CTF) determines the upper limit of data quality</li>
<li>Deep learning-based particle picking (Topaz/crYOLO) significantly improves performance on complex samples</li>
<li>2D classification screens particles; 3D classification handles conformational heterogeneity</li>
<li>Three key refinement steps: solvent mask + Bayesian polishing + CTF refinement</li>
<li>FSC = 0.143 is the gold standard for resolution</li>
</ul>
<p>Future articles will cover structure visualization (PyMOL/ChimeraX) and atomic modeling (Coot/phenix), turning density maps into atomic models.</p>`,Hb=`<h1>冷冻电镜单颗粒分析：数据处理全流程</h1>
<p>冷冻电镜单颗粒分析（Single Particle Analysis, SPA）是解析生物大分子近原子分辨率结构的核心技术。本文梳理从显微照片到原子模型的完整流程，介绍各步骤的目的、常用软件与关键参数。</p>
<h2>1. 流程总览</h2>
<pre><code>样品制备（冷冻制样）
   ↓
数据采集（Titan Krios / Glacios 等）
   ↓
① 预处理：运动校正 + CTF 估计
   ↓
② 颗粒挑选（Particle Picking）
   ↓
③ 2D 分类（筛除坏颗粒）
   ↓
④ 初始模型（Ab-initio）
   ↓
⑤ 3D 分类（异质性分析）
   ↓
⑥ 3D 精修（Refinement）
   ↓
⑦ 分辨率评估（FSC 曲线）
   ↓
⑧ 原子模型搭建与精修
</code></pre>
<h2>2. 样品制备与数据采集</h2>
<h3>2.1 冷冻制样</h3>
<ul>
<li>将蛋白样品（浓度通常 0.5–5 mg/mL）施加到载网（holey carbon grid）</li>
<li>用 Vitrobot / Leica GP 等设备快速浸入液乙烷，形成<strong>玻璃态冰</strong></li>
<li>理想冰厚：略薄于颗粒直径；冰太厚信噪比差，太薄颗粒易变形</li>
</ul>
<h3>2.2 数据采集</h3>
<ul>
<li>300 kV 电镜（Titan Krios）+ 直接电子探测器（Gatan K3 / Falcon 4）</li>
<li><strong>超分辨率模式</strong> + 剂量分配（total dose 约 40–60 e⁻/Å²）</li>
<li>每张照片（movie）通常 30–50 帧，用于运动校正</li>
<li>目标：10,000–20,000 张照片（高分辨率结构通常需要百万级颗粒）</li>
</ul>
<h2>3. 预处理（Preprocessing）</h2>
<h3>3.1 运动校正（Motion Correction）</h3>
<p>电子束引起的样品漂移会模糊图像，需要按帧校正：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RELION 4 中的 motion correction</span>
relion_run_motioncorr \\
  --i Micrographs/ \\
  --o MotionCorr/ \\
  --j 8 \\
  --pipeline_control MotionCorr/job001/
</code></pre>
<p>常用软件：<strong>MotionCor2</strong>（经典）、<strong>RELION 内置</strong>、<strong>cryoSPARC Patch Motion Correction</strong>。</p>
<h3>3.2 CTF 估计</h3>
<p>CTF（Contrast Transfer Function，衬度传递函数）描述了电子显微镜成像的相位翻转与欠焦效应：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># CTFFIND4（经典命令行）</span>
ctffind --input_micrograph MotionCorr/mic1.mrc \\
        --output_diag mic1_diag.mrc \\
        --defU 20000 --defV 20000

<span class="hljs-comment"># Gctf（更快的 GPU 版本）</span>
Gctf --input <span class="hljs-string">'MotionCorr/*.mrc'</span> --apix 1.0
</code></pre>
<p>常用软件：<strong>CTFFIND4</strong>、<strong>Gctf</strong>、<strong>cryoSPARC Patch CTF</strong>。</p>
<p><strong>关键判据</strong>：分辨率截断（Resolution limit）通常取 3–4 Å；剔除 astigmatism 过大或漂移严重的照片。</p>
<h2>4. 颗粒挑选（Particle Picking）</h2>
<h3>4.1 传统方法</h3>
<ul>
<li><strong>模板匹配</strong>（Template-based）：用 2D 类平均做模板，相关性搜索</li>
<li>软件：RELION 的 auto-picking、cryoSPARC 的 Template Picker</li>
</ul>
<h3>4.2 深度学习方法（当前主流）</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Topaz：基于 CNN 的颗粒挑选</span>
topaz train --model-path model.pkl --epochs 20 \\
  --num-workers 4 --train-data particles.toml

topaz extract --model-path model.pkl \\
  --micrographs micrographs.toml \\
  --output-prefix picked --radius 100
</code></pre>
<p>主流工具：<strong>Topaz</strong>、<strong>crYOLO</strong>（速度极快，支持实时挑选）、<strong>cryoSPARC Blob Picker</strong>。</p>
<blockquote>
<p>现代流程常先用 Blob/模板快速挑一批做 2D 分类，再用干净类平均训练 Topaz/crYOLO 模型，显著提升复杂样品（膜蛋白、小蛋白）的挑选质量。</p>
</blockquote>
<h2>5. 2D 分类：筛选颗粒</h2>
<p>2D 分类将颗粒按投影方向聚类，同时<strong>剔除冰污染、错误挑选、变性颗粒</strong>：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RELION 2D classification</span>
relion_refine --o Class2D/ \\
  --particle_images particles.star \\
  --ctf --K 100 --iter 25 --tau2_fudge 4 \\
  --flatten_solvent --zero_mask --pool 30 \\
  --per_particle_ctf --j 8
</code></pre>
<p><strong>质量判据</strong>：</p>
<ul>
<li>类平均应显示清晰的二级结构（α 螺旋、β 片层）</li>
<li>保留占总量 50–80% 的"好"颗粒</li>
<li>好的类平均是后续 3D 重构质量的前提</li>
</ul>
<h2>6. 3D 重构：初始模型与分类</h2>
<h3>6.1 初始模型（Ab-initio）</h3>
<p>无需参考模型，从随机起始重构：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># cryoSPARC</span>
cryosparc_abinitio ... --num_classes 3
</code></pre>
<ul>
<li>常用：<strong>cryoSPARC Ab-initio</strong>（速度快、稳健）、RELION 3D initial model</li>
<li>生成 2–4 个类，选择结构特征最清晰的类作为参考</li>
</ul>
<h3>6.2 3D 分类：处理构象异质性</h3>
<p>生物样品常存在多个构象（开放/闭合、配体结合/未结合）：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RELION 3D classification</span>
relion_refine --o Class3D/ \\
  --ref initial.mrc --K 4 --iter 40 \\
  --tau2_fudge 8 --particle_diameter 200
</code></pre>
<p><strong>策略</strong>：先粗分（K=4–6）看构象分布，再对感兴趣构象细分（K=2–3），最终每个 3D 类对应一个构象。</p>
<h2>7. 3D 精修（Refinement）</h2>
<h3>7.1 常规精修</h3>
<pre><code class="hljs language-bash">relion_refine --o Refine3D/ \\
  --ref class3d.mrc --particle_diameter 200 \\
  --auto_mask --solvent_mask \\
  --ctf --per_particle_ctf --j 8
</code></pre>
<p>关键参数：</p>
<ul>
<li><strong>C1 对称性</strong>起步，对称分子可指定 C/D/O 对称（如 C3、D7）</li>
<li><strong>溶剂掩膜</strong>（solvent mask）必须使用，否则 FSC 虚高</li>
<li>局部精修（local refinement）处理柔性区域</li>
</ul>
<h3>7.2 贝叶斯抛光（Bayesian Polishing）</h3>
<p>对 movie 帧级颗粒做轨迹拟合与加权，提升高频信息：</p>
<pre><code class="hljs language-bash">relion_polish --i Refine3D/run_data.star \\
  --o Polish/ \\
  --model_type 2 --nr_iter 5
</code></pre>
<h3>7.3 各向异性校正（CTF refinement）</h3>
<p>精修每个颗粒的欠焦、像散、束倾斜：</p>
<pre><code class="hljs language-bash">relion_ctf_refine --i Polish/particles.star \\
  --o CtfRefine/ \\
  --fit_defocus --fit_magnification --fit_tilt
</code></pre>
<h2>8. 分辨率评估：FSC 曲线</h2>
<p><strong>FSC（Fourier Shell Correlation）</strong>：把颗粒随机分成两个半集（half-maps），分别重构后比较两图的相关系数：</p>
<ul>
<li><strong>FSC = 0.143 准则</strong>：金标准分辨率判定（半图校正后）</li>
<li><strong>FSC = 0.5</strong>：保守估计（未校正）</li>
</ul>
<pre><code class="hljs language-bash">relion_postprocess --i Refine3D/run_half1_class001.mrc \\
  --mask mask.mrc --autob_masking \\
  --angpix 1.0
</code></pre>
<p><strong>质量检查表</strong>：</p>
<table>
<thead>
<tr>
<th>指标</th>
<th>良好标准</th>
</tr>
</thead>
<tbody>
<tr>
<td>全局分辨率（FSC 0.143）</td>
<td>≤ 3.5 Å（高分辨）</td>
</tr>
<tr>
<td>颗粒数</td>
<td>≥ 10 万（3 Å 级通常需数十万）</td>
</tr>
<tr>
<td>方向分布</td>
<td>均匀（无 Preferred orientation）</td>
</tr>
<tr>
<td>密度连续性</td>
<td>侧链密度清晰可见（&#x3C; 3 Å）</td>
</tr>
</tbody>
</table>
<h2>9. 常用软件栈对比</h2>
<table>
<thead>
<tr>
<th>环节</th>
<th>RELION（命令行）</th>
<th>cryoSPARC（GUI）</th>
</tr>
</thead>
<tbody>
<tr>
<td>运动校正</td>
<td>内置 / MotionCor2</td>
<td>Patch Motion</td>
</tr>
<tr>
<td>CTF</td>
<td>CTFFIND4 / Gctf</td>
<td>Patch CTF</td>
</tr>
<tr>
<td>挑选</td>
<td>模板 / Topaz 集成</td>
<td>Blob / Template / Topaz</td>
</tr>
<tr>
<td>2D/3D</td>
<td>稳定可靠</td>
<td>快、交互好</td>
</tr>
<tr>
<td>非均匀精修</td>
<td>—</td>
<td><strong>NU-refinement</strong>（柔性大分子利器）</td>
</tr>
</tbody>
</table>
<blockquote>
<p><strong>工作流建议</strong>：cryoSPARC 做预处理与初始重构，RELION 做高精度精修与抛光；或全流程 cryoSPARC（对新手更友好）。</p>
</blockquote>
<h2>10. 小结</h2>
<ul>
<li>预处理（运动校正 + CTF）决定数据质量上限</li>
<li>深度学习颗粒挑选（Topaz/crYOLO）显著提升复杂样品表现</li>
<li>2D 分类筛颗粒，3D 分类处理构象异质性</li>
<li>精修三件套：溶剂掩膜 + 贝叶斯抛光 + CTF refinement</li>
<li>FSC = 0.143 是分辨率金标准</li>
</ul>
<p>后续文章将介绍结构可视化（PyMOL/ChimeraX）与原子建模（Coot/phenix），把密度图变成原子模型。</p>`,Gb=`<h1>Atomic Modeling and Refinement</h1>
<p>After obtaining a high-quality density map, the next step is to build and refine an atomic model, ultimately yielding publishable, analyzable coordinate files. This article covers the complete workflow from automatic modeling to quality assessment.</p>
<h2>1. Workflow Overview</h2>
<pre><code>Density map (.mrc)
   ↓
① Initial model (automatic modeling / homologous model / AF prediction)
   ↓
② Model placement into density (docking / rigid body)
   ↓
③ Iterative refinement: Coot manual correction ↔ phenix automatic refinement
   ↓
④ Quality assessment (MolProbity / EMRinger / Q-score)
   ↓
⑤ Validation and publication (PDB deposition)
</code></pre>
<h2>2. Initial Model Sources</h2>
<h3>2.1 Automatic Modeling (Mainstream in the Deep Learning Era)</h3>
<ul>
<li><strong>ModelAngelo</strong>: Automatic modeling based on the AlphaFold architecture, building models directly from density maps</li>
<li><strong>DeepTracer</strong>: Fast fully automatic modeling (GPU)</li>
<li><strong>phenix.map_to_model</strong>: Classic automatic backbone/side-chain building</li>
</ul>
<pre><code class="hljs language-bash"><span class="hljs-comment"># ModelAngelo usage example</span>
model_angelo --halfmaps half1.mrc half2.mrc \\
  --output-dir model_angelo_out \\
  --num-workers 8

<span class="hljs-comment"># phenix.map_to_model</span>
phenix.map_to_model map.mrc resolution=3.0 \\
  output_model=auto_model.pdb
</code></pre>
<h3>2.2 Homology Modeling and Molecular Replacement</h3>
<ul>
<li>If a homologous structure exists: use <strong>phenix.dock_in_map</strong> to place it into the density</li>
<li>AlphaFold-predicted models: directly used as initial models (commonly used in cryo-EM)</li>
<li>At low resolution (> 5 Å), place rigid bodies first, then gradually refine</li>
</ul>
<h2>3. Coot: Manual Modeling and Correction</h2>
<p>Coot is the industry standard for manual model correction (GUI tool).</p>
<h3>3.1 Common Operations</h3>
<pre><code>File > Open Coordinates...   # Load model
File > Open Map...           # Load density map (FSC correlation calculated automatically)

# Basic corrections
1. Backbone adjustment: drag Cα (Real Space Refine Zone)
2. Side-chain orientation: Rotamers
3. Missing residues: Add Terminal Residue / Build
4. Delete extra: Delete Residue Range
5. Mutate residues: Mutate &#x26; Autofit
6. Ligand building: Ligand Builder / restraints
</code></pre>
<h3>3.2 Key Concept: Real Space Refine</h3>
<pre><code>Calculate > Real Space Refine Zone
</code></pre>
<ul>
<li>Select a density region, and Coot optimizes the geometry and density fit of that region</li>
<li>Repeat: fix → refine → check → fix again</li>
<li>After each refinement round, return to Coot to check problematic regions (Ramachandran outliers, clashes)</li>
</ul>
<h3>3.3 Ligand Modeling</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Generate ligand restraints</span>
phenix.elbow ligand.smiles    <span class="hljs-comment"># Generate chemical restraints from SMILES</span>
<span class="hljs-comment"># Or</span>
phenix.ligand_idealization ligand.cif
</code></pre>
<p>Use <code>Ligand > Ligand Builder</code> in Coot for manual building, then Real Space Refine.</p>
<h2>4. phenix Automatic Refinement</h2>
<h3>4.1 Basic Refinement</h3>
<pre><code class="hljs language-bash">phenix.real_space_refine model.pdb \\
  map.mrc \\
  resolution=3.0 \\
  output.prefix=refined \\
  <span class="hljs-built_in">nproc</span>=8
</code></pre>
<h3>4.2 Refinement Strategies</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Multiple rounds of iteration: global first, then local</span>
phenix.real_space_refine model.pdb map.mrc \\
  resolution=3.0 \\
  strategy=individual_sites+individual_adp \\
  restraints=ligand.cif \\
  <span class="hljs-built_in">nproc</span>=8
</code></pre>
<p>Commonly used strategy combinations:</p>
<ul>
<li><code>rigid_body</code>: rigid-body refinement first at low resolution</li>
<li><code>individual_sites</code>: per-atom coordinate refinement (high resolution)</li>
<li><code>individual_adp</code>: B-factor refinement</li>
<li><code>torsion</code>: torsion-angle refinement (commonly used for proteins)</li>
<li><code>local_grid</code>: local grid search (flexible regions)</li>
</ul>
<h3>4.3 Automatic Correction Loop (NCS and Water Molecules)</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Find water molecules (at high resolution)</span>
phenix.find_peaks_holes map.mrc map.mrc resolution=3.0

<span class="hljs-comment"># NCS restraints (significantly improves quality when symmetric subunits are present)</span>
phenix.real_space_refine ... strategy=individual_sites+individual_adp ncs_search.enabled=True
</code></pre>
<h2>5. Model Quality Assessment</h2>
<h3>5.1 Geometric Quality: MolProbity</h3>
<pre><code class="hljs language-bash">phenix.molprobity refined.pdb
</code></pre>
<p>Core metrics:</p>
<table>
<thead>
<tr>
<th>Metric</th>
<th>Good standard (3 Å)</th>
<th>Meaning</th>
</tr>
</thead>
<tbody>
<tr>
<td>Ramachandran favored</td>
<td>> 95%</td>
<td>Backbone dihedral angle validity</td>
</tr>
<tr>
<td>Ramachandran outliers</td>
<td>&#x3C; 0.5%</td>
<td>Outlier backbone</td>
</tr>
<tr>
<td>Rotamer outliers</td>
<td>&#x3C; 1%</td>
<td>Abnormal side-chain conformations</td>
</tr>
<tr>
<td>Clashscore</td>
<td>&#x3C; 5</td>
<td>Atomic collisions</td>
</tr>
<tr>
<td>MolProbity score</td>
<td>&#x3C; 2</td>
<td>Composite geometry score</td>
</tr>
</tbody>
</table>
<h3>5.2 Density Fit Quality</h3>
<ul>
<li><strong>EMRinger</strong>: assesses side-chain-to-density fit (> 2 indicates good)</li>
<li><strong>Q-score</strong>: assesses local density quality (1.0 = perfect, 0.5 = acceptable)</li>
<li><strong>CC_mask</strong>: model-density correlation coefficient (> 0.8 good)</li>
<li><strong>map-model FSC</strong>: FSC comparison between the model and the two half-maps</li>
</ul>
<pre><code class="hljs language-bash">phenix.map_model_cc refined.pdb map.mrc
phenix.em_ringer refined.pdb map.mrc
</code></pre>
<h3>5.3 Overfitting Detection</h3>
<p><strong>Free component (CC-free)</strong>: correlation computed independently from the half-map used in refinement. Must be monitored during refinement:</p>
<pre><code>Before refinement: CC_work ≈ CC_free
After refinement: CC_work significantly > CC_free → sign of overfitting
</code></pre>
<p><strong>Standard practice</strong>: use only one half-map for refinement (or full map + independent half-map validation), and use the other half-map to compute CC_free.</p>
<h2>6. Complete Refinement Workflow Example</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. Automatic modeling</span>
model_angelo --halfmaps half1.mrc half2.mrc --output-dir auto

<span class="hljs-comment"># 2. Initial refinement</span>
phenix.real_space_refine auto/model.pdb map.mrc \\
  resolution=3.1 strategy=rigid_body+individual_sites \\
  output.prefix=round1

<span class="hljs-comment"># 3. Coot manual correction (interactive)</span>
coot round1_refined.pdb map.mrc

<span class="hljs-comment"># 4. Second refinement round (add B-factors)</span>
phenix.real_space_refine coot_out.pdb map.mrc \\
  resolution=3.1 strategy=individual_sites+individual_adp \\
  output.prefix=round2

<span class="hljs-comment"># 5. Quality assessment</span>
phenix.molprobity round2_refined.pdb
phenix.map_model_cc round2_refined.pdb map.mrc
phenix.em_ringer round2_refined.pdb map.mrc

<span class="hljs-comment"># 6. Iterate until all metrics pass</span>
</code></pre>
<h2>7. PDB Deposition</h2>
<p>Before publication, prepare:</p>
<ol>
<li>Coordinate file (.pdb / .cif)</li>
<li>Density maps (half-maps + full map, .mrc)</li>
<li>Validation reports (MolProbity, EMRinger, etc.)</li>
<li>FSC curves and mask information</li>
</ol>
<p>Submit via <a href="https://deposit.wwpdb.org/">wwPDB deposition</a> to obtain PDB ID and EMDB ID.</p>
<h2>8. Common Problems</h2>
<table>
<thead>
<tr>
<th>Problem</th>
<th>Cause</th>
<th>Solution</th>
</tr>
</thead>
<tbody>
<tr>
<td>Blurred side-chain density</td>
<td>Insufficient resolution / incorrect B-factors</td>
<td>Increase resolution, check orientation distribution</td>
</tr>
<tr>
<td>Many Ramachandran outliers</td>
<td>Incorrect backbone building</td>
<td>Rebuild the region in Coot based on density</td>
</tr>
<tr>
<td>Overfitting</td>
<td>Excessive refinement</td>
<td>Monitor with CC_free, reduce iterations</td>
</tr>
<tr>
<td>Poor ligand fit</td>
<td>Incomplete restraints</td>
<td>Regenerate with phenix.elbow</td>
</tr>
<tr>
<td>Poor local region quality</td>
<td>Flexibility/conformational heterogeneity</td>
<td>Separate conformations with 3D classification, local refinement</td>
</tr>
</tbody>
</table>
<h2>9. Summary</h2>
<ul>
<li>Three sources for model building: automatic modeling (ModelAngelo/DeepTracer), homology/AF prediction, de novo building</li>
<li>Alternating iterations of Coot manual correction and phenix automatic refinement is the standard workflow</li>
<li>Triple quality check: geometry (MolProbity) + fit (EMRinger/Q-score) + overfitting (CC_free)</li>
<li>Sub-3 Å maps should aim for side-chain-level confidence; ligands and modifications require separate validation</li>
</ul>
<p>This completes the structural biology pipeline: data processing → technical review → visualization → atomic modeling, forming a complete loop from experiment to model.</p>`,Wb=`<h1>原子建模与精修</h1>
<p>拿到高质量密度图后，下一步是搭建原子模型并精修，最终获得可发布、可分析的坐标文件。本文覆盖从自动建模到质量评估的完整工作流。</p>
<h2>1. 工作流总览</h2>
<pre><code>密度图（.mrc）
   ↓
① 初始模型（自动建模 / 同源模型 / AF 预测）
   ↓
② 模型放入密度（docking / rigid body）
   ↓
③ 迭代精修：Coot 手动修正 ↔ phenix 自动精修
   ↓
④ 质量评估（MolProbity / EMRinger / Q-score）
   ↓
⑤ 验证与发布（PDB deposition）
</code></pre>
<h2>2. 初始模型来源</h2>
<h3>2.1 自动建模（深度学习时代的主流）</h3>
<ul>
<li><strong>ModelAngelo</strong>：基于 AlphaFold 架构的自动建模，直接从密度图搭模型</li>
<li><strong>DeepTracer</strong>：快速全自动建模（GPU）</li>
<li><strong>phenix.map_to_model</strong>：经典自动建主链/侧链</li>
</ul>
<pre><code class="hljs language-bash"><span class="hljs-comment"># ModelAngelo 使用示例</span>
model_angelo --halfmaps half1.mrc half2.mrc \\
  --output-dir model_angelo_out \\
  --num-workers 8

<span class="hljs-comment"># phenix.map_to_model</span>
phenix.map_to_model map.mrc resolution=3.0 \\
  output_model=auto_model.pdb
</code></pre>
<h3>2.2 同源建模与分子替换</h3>
<ul>
<li>有同源结构：用 <strong>phenix.dock_in_map</strong> 放入密度</li>
<li>AlphaFold 预测模型：直接作为初始模型（cryo-EM 中常用）</li>
<li>低分辨率（> 5 Å）时先放刚性体，再逐步细化</li>
</ul>
<h2>3. Coot：手动建模与修正</h2>
<p>Coot 是手动模型修正的行业标准（GUI 工具）。</p>
<h3>3.1 常用操作</h3>
<pre><code>File > Open Coordinates...   # 载入模型
File > Open Map...           # 载入密度图（自动算 FSC 相关）

# 基本修正
1. 主链调整：拖拽 Cα（Real Space Refine 区域）
2. 侧链方向：Rotamers（旋转异构体）
3. 缺失残基：Add Terminal Residue / Build
4. 删除多余：Delete Residue Range
5. 突变残基：Mutate &#x26; Autofit
6. 配体搭建：Ligand Builder / restraints
</code></pre>
<h3>3.2 关键理念：Real Space Refine</h3>
<pre><code>Calculate > Real Space Refine Zone
</code></pre>
<ul>
<li>选中密度区域，Coot 优化该区域几何与密度拟合</li>
<li>重复：修正 → 精修 → 检查 → 再修正</li>
<li>每轮精修后回到 Coot 检查不良区域（Ramachandran 离群点、clash）</li>
</ul>
<h3>3.3 配体建模</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 生成配体 restraints</span>
phenix.elbow ligand.smiles    <span class="hljs-comment"># 从 SMILES 生成化学约束</span>
<span class="hljs-comment"># 或</span>
phenix.ligand_idealization ligand.cif
</code></pre>
<p>Coot 中 <code>Ligand > Ligand Builder</code> 手动搭建，再 Real Space Refine。</p>
<h2>4. phenix 自动精修</h2>
<h3>4.1 基本精修</h3>
<pre><code class="hljs language-bash">phenix.real_space_refine model.pdb \\
  map.mrc \\
  resolution=3.0 \\
  output.prefix=refined \\
  <span class="hljs-built_in">nproc</span>=8
</code></pre>
<h3>4.2 精修策略</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 多轮迭代：先整体后局部</span>
phenix.real_space_refine model.pdb map.mrc \\
  resolution=3.0 \\
  strategy=individual_sites+individual_adp \\
  restraints=ligand.cif \\
  <span class="hljs-built_in">nproc</span>=8
</code></pre>
<p>常用 strategy 组合：</p>
<ul>
<li><code>rigid_body</code>：低分辨率先做刚性体</li>
<li><code>individual_sites</code>：逐原子坐标精修（高分辨率）</li>
<li><code>individual_adp</code>：B 因子精修</li>
<li><code>torsion</code>：扭转角精修（蛋白质常用）</li>
<li><code>local_grid</code>：局部网格搜索（柔性区域）</li>
</ul>
<h3>4.3 自动修正循环（NCS 与水分子）</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 找水分子（高分辨时）</span>
phenix.find_peaks_holes map.mrc map.mrc resolution=3.0

<span class="hljs-comment"># NCS 约束（有对称亚基时显著提升质量）</span>
phenix.real_space_refine ... strategy=individual_sites+individual_adp ncs_search.enabled=True
</code></pre>
<h2>5. 模型质量评估</h2>
<h3>5.1 几何质量：MolProbity</h3>
<pre><code class="hljs language-bash">phenix.molprobity refined.pdb
</code></pre>
<p>核心指标：</p>
<table>
<thead>
<tr>
<th>指标</th>
<th>良好标准（3 Å）</th>
<th>含义</th>
</tr>
</thead>
<tbody>
<tr>
<td>Ramachandran  favored</td>
<td>> 95%</td>
<td>主链二面角合理性</td>
</tr>
<tr>
<td>Ramachandran outliers</td>
<td>&#x3C; 0.5%</td>
<td>离群主链</td>
</tr>
<tr>
<td>Rotamer outliers</td>
<td>&#x3C; 1%</td>
<td>侧链构象异常</td>
</tr>
<tr>
<td>Clashscore</td>
<td>&#x3C; 5</td>
<td>原子碰撞</td>
</tr>
<tr>
<td>MolProbity score</td>
<td>&#x3C; 2</td>
<td>综合几何分数</td>
</tr>
</tbody>
</table>
<h3>5.2 密度拟合质量</h3>
<ul>
<li><strong>EMRinger</strong>：评估侧链与密度的拟合（> 2 表示良好）</li>
<li><strong>Q-score</strong>：评估局部密度质量（1.0 = 完美，0.5 = 可接受）</li>
<li><strong>CC_mask</strong>：模型-密度相关系数（> 0.8 良好）</li>
<li><strong>map-model FSC</strong>：模型与两张半图的 FSC 对比</li>
</ul>
<pre><code class="hljs language-bash">phenix.map_model_cc refined.pdb map.mrc
phenix.em_ringer refined.pdb map.mrc
</code></pre>
<h3>5.3 过拟合检测</h3>
<p><strong>Free component（CC-free）</strong>：与精修中使用的半图独立计算的相关性。精修时必须监控：</p>
<pre><code>精修前 CC_work ≈ CC_free
精修后 CC_work 明显 > CC_free → 过拟合信号
</code></pre>
<p><strong>规范做法</strong>：精修只用一张半图（或全图 + 独立半图验证），用另一张半图计算 CC_free。</p>
<h2>6. 完整精修流程示例</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. 自动建模</span>
model_angelo --halfmaps half1.mrc half2.mrc --output-dir auto

<span class="hljs-comment"># 2. 初始精修</span>
phenix.real_space_refine auto/model.pdb map.mrc \\
  resolution=3.1 strategy=rigid_body+individual_sites \\
  output.prefix=round1

<span class="hljs-comment"># 3. Coot 手动修正（交互）</span>
coot round1_refined.pdb map.mrc

<span class="hljs-comment"># 4. 第二轮精修（加 B 因子）</span>
phenix.real_space_refine coot_out.pdb map.mrc \\
  resolution=3.1 strategy=individual_sites+individual_adp \\
  output.prefix=round2

<span class="hljs-comment"># 5. 质量评估</span>
phenix.molprobity round2_refined.pdb
phenix.map_model_cc round2_refined.pdb map.mrc
phenix.em_ringer round2_refined.pdb map.mrc

<span class="hljs-comment"># 6. 迭代直到所有指标达标</span>
</code></pre>
<h2>7. 提交 PDB</h2>
<p>发布前需准备：</p>
<ol>
<li>坐标文件（.pdb / .cif）</li>
<li>密度图（half-maps + full map，.mrc）</li>
<li>验证报告（MolProbity、EMRinger 等）</li>
<li>FSC 曲线与 mask 信息</li>
</ol>
<p>通过 <a href="https://deposit.wwpdb.org/">wwPDB deposition</a> 提交，获取 PDB ID 与 EMDB ID。</p>
<h2>8. 常见问题</h2>
<table>
<thead>
<tr>
<th>问题</th>
<th>原因</th>
<th>解决</th>
</tr>
</thead>
<tbody>
<tr>
<td>侧链密度模糊</td>
<td>分辨率不足 / B 因子错误</td>
<td>提高分辨率、检查方向分布</td>
</tr>
<tr>
<td>Ramachandran 离群多</td>
<td>主链建错</td>
<td>Coot 中按密度重建该区域</td>
</tr>
<tr>
<td>过拟合</td>
<td>精修过度</td>
<td>用 CC_free 监控、减少迭代</td>
</tr>
<tr>
<td>配体拟合差</td>
<td>restraints 不完整</td>
<td>phenix.elbow 重新生成</td>
</tr>
<tr>
<td>局部区域质量差</td>
<td>柔性/构象异质性</td>
<td>3D 分类分离构象、局部精修</td>
</tr>
</tbody>
</table>
<h2>9. 小结</h2>
<ul>
<li>建模三来源：自动建模（ModelAngelo/DeepTracer）、同源/AF 预测、从头搭建</li>
<li>Coot 手动修正与 phenix 自动精修<strong>交替迭代</strong>是标准节奏</li>
<li>质量三重检查：几何（MolProbity）+ 拟合（EMRinger/Q-score）+ 过拟合（CC_free）</li>
<li>亚 3 Å 应追求侧链级可信度，配体与修饰需单独验证</li>
</ul>
<p>至此结构生物学方向完成：数据处理流程 → 技术综述 → 可视化 → 原子建模，构成从实验到模型的完整闭环。</p>`,Kb=`<h1>Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice</h1>
<p>Structure visualization is an everyday tool in structural biology. Using two major software packages, PyMOL and UCSF ChimeraX, this article covers everything from PDB data retrieval to publication-quality rendering, encompassing all the core operations needed in daily research.</p>
<h2>1. Data Source: The PDB Database</h2>
<h3>1.1 Search and Download</h3>
<p><a href="https://www.rcsb.org/">wwPDB / RCSB PDB</a> is the global repository for structural data:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Command-line download (by PDB ID)</span>
wget https://files.rcsb.org/download/1CRN.pdb
wget https://files.rcsb.org/download/7A0A.pdb   <span class="hljs-comment"># Cryo-EM structure example</span>
</code></pre>
<h3>1.2 Core Contents of a PDB File</h3>
<pre><code>HEADER    PLANT SEED PROTEIN             08-JAN-81   1CRN
TITLE     CRAMBIN
ATOM      1  N   THR A   1      11.106  16.700  17.083  1.00 20.14           N
ATOM      2  CA  THR A   1      11.579  17.437  17.166  1.00 20.55           C
...
HELIX    1   1 THR A    1  SER A    8  1                                  8
SHEET    1   A 2 ILE A   9  ARG A  12  0
</code></pre>
<ul>
<li><code>ATOM</code> lines: atomic coordinates (residue, chain, x/y/z, B-factor, element)</li>
<li><code>HELIX</code> / <code>SHEET</code>: secondary structure annotations</li>
<li><code>CONECT</code>: covalent connections</li>
<li>The newer <code>mmCIF</code> format (<code>.cif</code>) has become mainstream and carries more complete information</li>
</ul>
<h2>2. PyMOL Basics</h2>
<h3>2.1 Launch and Load</h3>
<pre><code class="hljs language-bash">pymol 1CRN.pdb          <span class="hljs-comment"># Load directly</span>
<span class="hljs-comment"># Or after launching</span>
<span class="hljs-comment"># File > Open</span>
</code></pre>
<h3>2.2 Basic Display Commands</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Load</span>
fetch 1crn               <span class="hljs-comment"># Fetch from the network</span>
load 1CRN.pdb

<span class="hljs-comment"># Display modes</span>
show cartoon             <span class="hljs-comment"># Cartoon (backbone trace)</span>
show sticks              <span class="hljs-comment"># Sticks (atomic details)</span>
show surface             <span class="hljs-comment"># Surface</span>
hide lines               <span class="hljs-comment"># Hide lines</span>

<span class="hljs-comment"># Coloring</span>
color cyan               <span class="hljs-comment"># Overall color</span>
spectrum count           <span class="hljs-comment"># Rainbow colors (by residue number)</span>
color red, resi <span class="hljs-number">1</span>-<span class="hljs-number">10</span>     <span class="hljs-comment"># Color a specific region</span>
</code></pre>
<h3>2.3 Selections and Operations</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Selection syntax</span>
select helix, ss h                   <span class="hljs-comment"># All helices</span>
select sheet, ss s                   <span class="hljs-comment"># All beta sheets</span>
select active_site, resi <span class="hljs-number">15</span>+<span class="hljs-number">20</span>+<span class="hljs-number">45</span>    <span class="hljs-comment"># Specific residues</span>
select ligand, resn HEM              <span class="hljs-comment"># Ligand (heme)</span>

<span class="hljs-comment"># Object operations</span>
show sticks, active_site
color yellow, active_site
zoom active_site
</code></pre>
<h3>2.4 Measurements</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Distance</span>
distance d1, /1CRN//A/<span class="hljs-number">15</span>/CA, /1CRN//A/<span class="hljs-number">20</span>/CA

<span class="hljs-comment"># Angle / dihedral</span>
angle a1, /1CRN//A/<span class="hljs-number">15</span>/CA, /1CRN//A/<span class="hljs-number">16</span>/CA, /1CRN//A/<span class="hljs-number">17</span>/CA
dihedral dh1, /1CRN//A/<span class="hljs-number">15</span>/CA, /1CRN//A/<span class="hljs-number">16</span>/CA, /1CRN//A/<span class="hljs-number">17</span>/CA, /1CRN//A/<span class="hljs-number">18</span>/CA
</code></pre>
<h3>2.5 Interaction Analysis</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Hydrogen bonds</span>
distance hbond, (resn HEM), (resn HIS)

<span class="hljs-comment"># Contact interface</span>
select interface, chain A within <span class="hljs-number">4.5</span> of chain B
</code></pre>
<h2>3. PyMOL Script Batch Processing</h2>
<p>Write commonly used operations as scripts for reproducible execution:</p>
<pre><code class="hljs language-python"><span class="hljs-comment"># render_view.py</span>
fetch 1crn
hide everything
show cartoon
color spectrum, ss
<span class="hljs-built_in">set</span> cartoon_transparency, <span class="hljs-number">0.2</span>

<span class="hljs-comment"># Zoom into the binding site and render</span>
zoom resi <span class="hljs-number">1</span>-<span class="hljs-number">20</span>
<span class="hljs-built_in">set</span> ray_shadows, <span class="hljs-number">1</span>
<span class="hljs-built_in">set</span> ray_opaque_background, <span class="hljs-number">0</span>
ray <span class="hljs-number">1200</span>, <span class="hljs-number">900</span>
png 1crn_view.png, dpi=<span class="hljs-number">300</span>
</code></pre>
<pre><code class="hljs language-bash">pymol -cq render_view.py    <span class="hljs-comment"># -c command-line mode, -q quiet</span>
</code></pre>
<h2>4. UCSF ChimeraX Basics</h2>
<p>ChimeraX has a friendly, modern interface (Qt GUI) and supports batch operations via the <code>open</code> command:</p>
<h3>4.1 Opening Structures</h3>
<pre><code class="hljs language-bash">chimerax 1CRN.pdb
<span class="hljs-comment"># Or via command</span>
open 1crn
open 7a0a
</code></pre>
<h3>4.2 Common Commands</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Display modes</span>
cartoon
stick
surface

<span class="hljs-comment"># Coloring</span>
color bychain          <span class="hljs-comment"># Color by chain</span>
color byhetero         <span class="hljs-comment"># Distinguish ligands/ions</span>
color byattribute bfactor   <span class="hljs-comment"># Color by B-factor</span>

<span class="hljs-comment"># Selections</span>
select :<span class="hljs-number">15</span>-<span class="hljs-number">20</span>          <span class="hljs-comment"># Residues 15-20</span>
select /A               <span class="hljs-comment"># Chain A</span>
select ligand
</code></pre>
<h3>4.3 Cryo-EM Density Map Display (ChimeraX's Strength)</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Open a density map (.mrc)</span>
<span class="hljs-built_in">open</span> <span class="hljs-built_in">map</span>.mrc

<span class="hljs-comment"># Adjust isosurface level</span>
volume level <span class="hljs-number">0.01</span>
volume <span class="hljs-comment">#2 level 0.05</span>

<span class="hljs-comment"># Superpose density map and model</span>
<span class="hljs-built_in">open</span> model.pdb
<span class="hljs-built_in">open</span> <span class="hljs-built_in">map</span>.mrc
volume <span class="hljs-comment">#2 level 0.02</span>
color <span class="hljs-comment">#1 cornflowerblue</span>
transparency <span class="hljs-comment">#2 30</span>

<span class="hljs-comment"># Local density inspection</span>
volume zone <span class="hljs-comment">#2 near :45-60 radius 5</span>
</code></pre>
<p><strong>ChimeraX's <code>volume zone</code> is the core tool for checking local density quality</strong>: display density around residues to assess whether side chains are discernible.</p>
<h2>5. Key Points for High-Quality Rendering</h2>
<h3>5.1 Color Principles</h3>
<ul>
<li>Cartoon coloring: by chain (<code>bychain</code>) or by domain</li>
<li>Key residues: use a single highlight color (yellow/red/blue), avoid overusing rainbows</li>
<li>Ligands: <code>byhetero</code> automatically assigns element-based colors (C gray, N blue, O red, S yellow)</li>
</ul>
<h3>5.2 Lighting and Materials</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># PyMOL</span>
<span class="hljs-built_in">set</span> ray_shadows, <span class="hljs-number">1</span>
<span class="hljs-built_in">set</span> specular, <span class="hljs-number">0.5</span>
<span class="hljs-built_in">set</span> ambient, <span class="hljs-number">0.3</span>

<span class="hljs-comment"># ChimeraX</span>
<span class="hljs-built_in">set</span> lightMode full
graphics silhouettes true
</code></pre>
<h3>5.3 Resolution Output</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># PyMOL publication-quality rendering</span>
ray 2400, 1800
png figure.png, dpi=600

<span class="hljs-comment"># ChimeraX</span>
save figure.png width 2400 height 1800 supersample 3
</code></pre>
<h2>6. Structure Comparison and Superposition</h2>
<pre><code class="hljs language-python"><span class="hljs-comment"># PyMOL: align two structures</span>
align model2, model1

<span class="hljs-comment"># ChimeraX: matchmaker</span>
matchmaker <span class="hljs-comment">#2 to #1</span>
</code></pre>
<p><strong>RMSD</strong> is the standard metric for measuring structural similarity (typically reported as Cα or all-atom RMSD).</p>
<h2>7. Interactive Viewing and Sharing</h2>
<ul>
<li><em><em>Mol</em> / NGL</em>*: web-based 3D viewing (built into PDBe, RCSB)</li>
<li><strong>PyMOL Web</strong>: export sessions as interactive HTML pages</li>
<li>Commonly used for scientific communication: PyMOL sessions (<code>.pse</code>), ChimeraX sessions (<code>.cxs</code>)</li>
</ul>
<h2>8. Summary</h2>
<ul>
<li>PDB is the data source; mmCIF is the new standard</li>
<li>PyMOL: powerful commands, scriptable batch processing (<code>pymol -cq script.py</code>)</li>
<li>ChimeraX: modern GUI + density map handling (<code>volume</code>, <code>volume zone</code>) is irreplaceable</li>
<li>Publication-quality rendering: hide the clutter → highlight the key → high-quality lighting → high-resolution export</li>
</ul>
<p>The next article will cover atomic modeling: how to build and refine atomic models from density maps.</p>`,Xb=`<h1>生物大分子结构可视化：PyMOL 与 ChimeraX 实战</h1>
<p>结构可视化是结构生物学的日常工具。本文以 PyMOL 与 UCSF ChimeraX 两大主流软件为例，从 PDB 数据获取到出版级渲染，覆盖日常科研所需的全部核心操作。</p>
<h2>1. 数据来源：PDB 数据库</h2>
<h3>1.1 检索与下载</h3>
<p><a href="https://www.rcsb.org/">wwPDB / RCSB PDB</a> 是全球结构数据仓库：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 命令行下载（按 PDB ID）</span>
wget https://files.rcsb.org/download/1CRN.pdb
wget https://files.rcsb.org/download/7A0A.pdb   <span class="hljs-comment"># 冷冻电镜结构示例</span>
</code></pre>
<h3>1.2 PDB 文件的核心内容</h3>
<pre><code>HEADER    PLANT SEED PROTEIN             08-JAN-81   1CRN
TITLE     CRAMBIN
ATOM      1  N   THR A   1      11.106  16.700  17.083  1.00 20.14           N
ATOM      2  CA  THR A   1      11.579  17.437  17.166  1.00 20.55           C
...
HELIX    1   1 THR A    1  SER A    8  1                                  8
SHEET    1   A 2 ILE A   9  ARG A  12  0
</code></pre>
<ul>
<li><code>ATOM</code> 行：原子坐标（残基、链、x/y/z、B 因子、元素）</li>
<li><code>HELIX</code> / <code>SHEET</code>：二级结构注释</li>
<li><code>CONECT</code>：共价连接</li>
<li>新格式 <code>mmCIF</code>（.cif）已成为主流，信息更完整</li>
</ul>
<h2>2. PyMOL 基础</h2>
<h3>2.1 启动与加载</h3>
<pre><code class="hljs language-bash">pymol 1CRN.pdb          <span class="hljs-comment"># 直接加载</span>
<span class="hljs-comment"># 或启动后</span>
<span class="hljs-comment"># File > Open</span>
</code></pre>
<h3>2.2 基础显示命令</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 载入</span>
fetch 1crn               <span class="hljs-comment"># 从网络获取</span>
load 1CRN.pdb

<span class="hljs-comment"># 显示模式</span>
show cartoon             <span class="hljs-comment"># 卡通（主链走向）</span>
show sticks              <span class="hljs-comment"># 棍状（原子细节）</span>
show surface             <span class="hljs-comment"># 表面</span>
hide lines               <span class="hljs-comment"># 隐藏线条</span>

<span class="hljs-comment"># 着色</span>
color cyan               <span class="hljs-comment"># 整体着色</span>
spectrum count           <span class="hljs-comment"># 彩虹色（按残基序号）</span>
color red, resi <span class="hljs-number">1</span>-<span class="hljs-number">10</span>     <span class="hljs-comment"># 指定区域着色</span>
</code></pre>
<h3>2.3 选择与操作</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 选择语法</span>
select helix, ss h                   <span class="hljs-comment"># 所有螺旋</span>
select sheet, ss s                   <span class="hljs-comment"># 所有片层</span>
select active_site, resi <span class="hljs-number">15</span>+<span class="hljs-number">20</span>+<span class="hljs-number">45</span>    <span class="hljs-comment"># 指定残基</span>
select ligand, resn HEM              <span class="hljs-comment"># 配体（血红素）</span>

<span class="hljs-comment"># 对象操作</span>
show sticks, active_site
color yellow, active_site
zoom active_site
</code></pre>
<h3>2.4 测量</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 距离</span>
distance d1, /1CRN//A/<span class="hljs-number">15</span>/CA, /1CRN//A/<span class="hljs-number">20</span>/CA

<span class="hljs-comment"># 角度 / 二面角</span>
angle a1, /1CRN//A/<span class="hljs-number">15</span>/CA, /1CRN//A/<span class="hljs-number">16</span>/CA, /1CRN//A/<span class="hljs-number">17</span>/CA
dihedral dh1, /1CRN//A/<span class="hljs-number">15</span>/CA, /1CRN//A/<span class="hljs-number">16</span>/CA, /1CRN//A/<span class="hljs-number">17</span>/CA, /1CRN//A/<span class="hljs-number">18</span>/CA
</code></pre>
<h3>2.5 相互作用分析</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 氢键</span>
distance hbond, (resn HEM), (resn HIS)

<span class="hljs-comment"># 接触面</span>
select interface, chain A within <span class="hljs-number">4.5</span> of chain B
</code></pre>
<h2>3. PyMOL 脚本批处理</h2>
<p>把常用操作写成脚本，可重复执行：</p>
<pre><code class="hljs language-python"><span class="hljs-comment"># render_view.py</span>
fetch 1crn
hide everything
show cartoon
color spectrum, ss
<span class="hljs-built_in">set</span> cartoon_transparency, <span class="hljs-number">0.2</span>

<span class="hljs-comment"># 放大结合位点并渲染</span>
zoom resi <span class="hljs-number">1</span>-<span class="hljs-number">20</span>
<span class="hljs-built_in">set</span> ray_shadows, <span class="hljs-number">1</span>
<span class="hljs-built_in">set</span> ray_opaque_background, <span class="hljs-number">0</span>
ray <span class="hljs-number">1200</span>, <span class="hljs-number">900</span>
png 1crn_view.png, dpi=<span class="hljs-number">300</span>
</code></pre>
<pre><code class="hljs language-bash">pymol -cq render_view.py    <span class="hljs-comment"># -c 命令行模式，-q 静默</span>
</code></pre>
<h2>4. UCSF ChimeraX 基础</h2>
<p>ChimeraX 界面友好、现代（Qt GUI），且支持 <code>open</code> 命令批量操作：</p>
<h3>4.1 打开结构</h3>
<pre><code class="hljs language-bash">chimerax 1CRN.pdb
<span class="hljs-comment"># 或命令</span>
open 1crn
open 7a0a
</code></pre>
<h3>4.2 常用命令</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 显示模式</span>
cartoon
stick
surface

<span class="hljs-comment"># 配色</span>
color bychain          <span class="hljs-comment"># 按链着色</span>
color byhetero         <span class="hljs-comment"># 配体/离子区分</span>
color byattribute bfactor   <span class="hljs-comment"># 按 B 因子</span>

<span class="hljs-comment"># 选择</span>
select :<span class="hljs-number">15</span>-<span class="hljs-number">20</span>          <span class="hljs-comment"># 残基 15-20</span>
select /A               <span class="hljs-comment"># A 链</span>
select ligand
</code></pre>
<h3>4.3 冷冻电镜密度图展示（ChimeraX 强项）</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 打开密度图（.mrc）</span>
<span class="hljs-built_in">open</span> <span class="hljs-built_in">map</span>.mrc

<span class="hljs-comment"># 调整等值面水平</span>
volume level <span class="hljs-number">0.01</span>
volume <span class="hljs-comment">#2 level 0.05</span>

<span class="hljs-comment"># 密度图与模型叠加</span>
<span class="hljs-built_in">open</span> model.pdb
<span class="hljs-built_in">open</span> <span class="hljs-built_in">map</span>.mrc
volume <span class="hljs-comment">#2 level 0.02</span>
color <span class="hljs-comment">#1 cornflowerblue</span>
transparency <span class="hljs-comment">#2 30</span>

<span class="hljs-comment"># 区域密度查看（局部）</span>
volume zone <span class="hljs-comment">#2 near :45-60 radius 5</span>
</code></pre>
<p><strong>ChimeraX 的 <code>volume zone</code> 是检查局部密度质量的核心工具</strong>：在残基周围显示密度，判断侧链是否可辨。</p>
<h2>5. 高质量渲染要点</h2>
<h3>5.1 色彩原则</h3>
<ul>
<li>卡通着色：按链（bychain）或按结构域</li>
<li>关键残基：统一高亮色（黄/红/蓝），避免彩虹滥用</li>
<li>配体：<code>byhetero</code> 自动区分元素色（C 灰、N 蓝、O 红、S 黄）</li>
</ul>
<h3>5.2 光照与材质</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># PyMOL</span>
<span class="hljs-built_in">set</span> ray_shadows, <span class="hljs-number">1</span>
<span class="hljs-built_in">set</span> specular, <span class="hljs-number">0.5</span>
<span class="hljs-built_in">set</span> ambient, <span class="hljs-number">0.3</span>

<span class="hljs-comment"># ChimeraX</span>
<span class="hljs-built_in">set</span> lightMode full
graphics silhouettes true
</code></pre>
<h3>5.3 分辨率输出</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># PyMOL 出版级渲染</span>
ray 2400, 1800
png figure.png, dpi=600

<span class="hljs-comment"># ChimeraX</span>
save figure.png width 2400 height 1800 supersample 3
</code></pre>
<h2>6. 结构比对与叠合</h2>
<pre><code class="hljs language-python"><span class="hljs-comment"># PyMOL：对齐两个结构</span>
align model2, model1

<span class="hljs-comment"># ChimeraX：matchmaker</span>
matchmaker <span class="hljs-comment">#2 to #1</span>
</code></pre>
<p><strong>RMSD</strong> 是衡量结构相似度的标准指标（通常报告 Cα 或全原子 RMSD）。</p>
<h2>7. 交互式查看与分享</h2>
<ul>
<li><em><em>Mol</em> / NGL</em>*：网页端 3D 查看（PDBe、RCSB 内置）</li>
<li><strong>PyMOL Web</strong>：把会话导出为 HTML 交互页面</li>
<li>科研交流常用：PyMOL session（.pse）、ChimeraX session（.cxs）</li>
</ul>
<h2>8. 小结</h2>
<ul>
<li>PDB 是数据源头，mmCIF 是新标准</li>
<li>PyMOL：命令强大、脚本化批处理（<code>pymol -cq script.py</code>）</li>
<li>ChimeraX：现代 GUI + 密度图处理（<code>volume</code>、<code>volume zone</code>）无可替代</li>
<li>出版级渲染：隐藏杂项 → 突出关键 → 高质量光线 → 高分辨率导出</li>
</ul>
<p>下一篇将介绍原子建模：如何从密度图搭建并精修原子模型。</p>`,Nt=()=>{const s=globalThis.__GBLOG_LOCALE__;return s||(typeof window<"u"?localStorage.getItem("locale"):null)||"zh-CN"},Yb=(s,n=".json")=>{const e=Nt()==="en-US";return s.endsWith("-en")?`${s}${n}`:e?`${s}-en${n}`:`${s}${n}`},Xe=Object.assign({"/data-branch/content/about-en.json":Hg,"/data-branch/content/about.json":Wg,"/data-branch/content/categories-en.json":Xg,"/data-branch/content/categories.json":Qg,"/data-branch/content/notes-en.json":Zg,"/data-branch/content/notes.json":nb,"/data-branch/content/posts-en.json":eb,"/data-branch/content/posts.json":lb,"/data-branch/content/projects-en.json":pb,"/data-branch/content/projects.json":rb,"/data-branch/content/resources-en.json":ub,"/data-branch/content/resources.json":db,"/data-branch/content/tags-en.json":jb,"/data-branch/content/tags.json":gb,"/data-branch/content/topics-en.json":_b,"/data-branch/content/topics.json":vb}),mc=Object.assign({"/data-branch/content/html/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en.html":wb,"/data-branch/content/html/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools.html":Cb,"/data-branch/content/html/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools-en.html":kb,"/data-branch/content/html/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools.html":xb,"/data-branch/content/html/notes/Programming/bash/bash-scripting/bash-scripting-en.html":Pb,"/data-branch/content/html/notes/Programming/bash/bash-scripting/bash-scripting.html":Eb,"/data-branch/content/html/notes/Programming/linux/linux-basics/linux-basics-en.html":Sb,"/data-branch/content/html/notes/Programming/linux/linux-basics/linux-basics.html":Tb,"/data-branch/content/html/notes/Programming/python/python-advanced/python-advanced-en.html":Ab,"/data-branch/content/html/notes/Programming/python/python-advanced/python-advanced.html":Rb,"/data-branch/content/html/notes/Programming/python/python-basics/python-basics-en.html":Lb,"/data-branch/content/html/notes/Programming/python/python-basics/python-basics.html":Mb,"/data-branch/content/html/notes/Programming/python/python-data/python-data-en.html":Db,"/data-branch/content/html/notes/Programming/python/python-data/python-data.html":Ob,"/data-branch/content/html/notes/Programming/r/r-basics/r-basics-en.html":Ib,"/data-branch/content/html/notes/Programming/r/r-basics/r-basics.html":Nb,"/data-branch/content/html/notes/Programming/r/r-ggplot2/r-ggplot2-en.html":Fb,"/data-branch/content/html/notes/Programming/r/r-ggplot2/r-ggplot2.html":$b,"/data-branch/content/html/notes/Programming/r/r-tidyverse/r-tidyverse-en.html":Bb,"/data-branch/content/html/notes/Programming/r/r-tidyverse/r-tidyverse.html":qb,"/data-branch/content/html/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en.html":zb,"/data-branch/content/html/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview.html":Vb,"/data-branch/content/html/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow-en.html":Ub,"/data-branch/content/html/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow.html":Hb,"/data-branch/content/html/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling-en.html":Gb,"/data-branch/content/html/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling.html":Wb,"/data-branch/content/html/notes/StructuralBiology/visualization/structure-visualization/structure-visualization-en.html":Kb,"/data-branch/content/html/notes/StructuralBiology/visualization/structure-visualization/structure-visualization.html":Xb}),jc=Object.assign({"/data-branch/cache/en/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en.md":()=>ys(()=>import("./protein-design-tools-en-kwRRVDCZ.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools-en.md":()=>ys(()=>import("./virtual-screening-tools-en-Ds1a_3n3.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/bash/bash-scripting/bash-scripting-en.md":()=>ys(()=>import("./bash-scripting-en-DWq6sXnJ.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/linux/linux-basics/linux-basics-en.md":()=>ys(()=>import("./linux-basics-en-DSb4s4PH.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/python/python-advanced/python-advanced-en.md":()=>ys(()=>import("./python-advanced-en-BQKRq7MH.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/python/python-basics/python-basics-en.md":()=>ys(()=>import("./python-basics-en-BvTX6dVV.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/python/python-data/python-data-en.md":()=>ys(()=>import("./python-data-en-DTE7E1vm.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/r/r-basics/r-basics-en.md":()=>ys(()=>import("./r-basics-en-Dzr6Fzbi.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/r/r-ggplot2/r-ggplot2-en.md":()=>ys(()=>import("./r-ggplot2-en-CvP2vjkk.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/r/r-tidyverse/r-tidyverse-en.md":()=>ys(()=>import("./r-tidyverse-en-CuLVIm9c.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en.md":()=>ys(()=>import("./cryoem-overview-en-DsP4yYlt.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow-en.md":()=>ys(()=>import("./cryoem-workflow-en-BpFgioRH.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling-en.md":()=>ys(()=>import("./atomic-modeling-en-BX1y6yZD.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/StructuralBiology/visualization/structure-visualization/structure-visualization-en.md":()=>ys(()=>import("./structure-visualization-en-CHXUkj4Z.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en.md":()=>ys(()=>import("./protein-design-tools-en-BdlMbrNY.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools.md":()=>ys(()=>import("./protein-design-tools-Jzs48-a0.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools-en.md":()=>ys(()=>import("./virtual-screening-tools-en-BEtVbCe5.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools.md":()=>ys(()=>import("./virtual-screening-tools-BdHWOXdu.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/bash/bash-scripting/bash-scripting-en.md":()=>ys(()=>import("./bash-scripting-en-CYovvjQL.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/bash/bash-scripting/bash-scripting.md":()=>ys(()=>import("./bash-scripting-DEkzGDwa.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/linux/linux-basics/linux-basics-en.md":()=>ys(()=>import("./linux-basics-en-BDL0HlS2.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/linux/linux-basics/linux-basics.md":()=>ys(()=>import("./linux-basics-UyN4mB_u.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-advanced/python-advanced-en.md":()=>ys(()=>import("./python-advanced-en-D3xdyM1A.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-advanced/python-advanced.md":()=>ys(()=>import("./python-advanced-iN1UORRl.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-basics/python-basics-en.md":()=>ys(()=>import("./python-basics-en-DOB3GCxE.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-basics/python-basics.md":()=>ys(()=>import("./python-basics-uh9Cz9FQ.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-data/python-data-en.md":()=>ys(()=>import("./python-data-en-BH1w0rP8.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-data/python-data.md":()=>ys(()=>import("./python-data-CmFq1YcA.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/r/r-basics/r-basics-en.md":()=>ys(()=>import("./r-basics-en-CPaztwoM.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/r/r-basics/r-basics.md":()=>ys(()=>import("./r-basics-CqBb8X8p.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/r/r-ggplot2/r-ggplot2.md":()=>ys(()=>import("./r-ggplot2-CK3q5J49.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/r/r-tidyverse/r-tidyverse-en.md":()=>ys(()=>import("./r-tidyverse-en-C13YYJSL.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/r/r-tidyverse/r-tidyverse.md":()=>ys(()=>import("./r-tidyverse-v7eXFCz-.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview.md":()=>ys(()=>import("./cryoem-overview-DT2Crx5Y.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow.md":()=>ys(()=>import("./cryoem-workflow-C2zs9nW6.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling.md":()=>ys(()=>import("./atomic-modeling-BGRUqNS4.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/StructuralBiology/visualization/structure-visualization/structure-visualization.md":()=>ys(()=>import("./structure-visualization-B3-STJgI.js"),[]).then(s=>s.default)}),Xa=s=>{const n=Yb(s,".json"),a=Object.keys(Xe).find(e=>e.includes(n));if(a)return Xe[a].default||{};if(console.error(`Failed to load JSON content: ${n}`),Nt()==="en-US"){const e=Object.keys(Xe).find(t=>t.includes(`${s}.json`));if(e)return Xe[e].default||{}}return{}},Qb=s=>{const a=Nt()==="en-US",e=s.replace(/\.md$/i,""),t=[];e.endsWith("-en")?(t.push(e),a&&t.push(e.replace(/-en$/,""))):(t.push(a?`${e}-en`:e),t.push(e));for(const l of t){const o=Object.keys(mc).find(p=>p.includes(`/content/html/${l}.html`));if(o)return mc[o]}return`<h1>Content Not Available</h1>
<p>The requested content could not be loaded.</p>`},Ft=()=>Xa("categories"),Jb=()=>Xa("posts"),Zb=()=>Xa("notes"),s_=()=>Xa("tags"),n_=()=>Xa("about"),a_=()=>Xa("resources"),e_=async s=>{const n=Nt()==="en-US"?"-en":"",a=`${s.replace(/\.md$/,"")}${n}.md`,e=Object.keys(jc).find(t=>t.endsWith(`/${a}`)||t.endsWith(a));if(!e)return null;try{return await jc[e]()}catch(t){return console.error(`Failed to load markdown source: ${a}`,t),null}},t_={class:"navigation-tree"},l_={class:"tree-label"},o_={class:"article-list article-list-root"},p_={key:0,class:"article-item"},c_={class:"directory-node"},r_={class:"tree-item tree-item--folder"},i_={key:0,class:"article-list sub-list files-level"},u_={key:1,class:"article-list sub-list"},h_={class:"directory-node"},d_={class:"tree-item tree-item--folder"},m_={key:0,class:"article-list sub-list files-level"},j_=Gs({__name:"NavigationTree",setup(s){const{locale:n}=nn(),a=Ka(),e=os([]),t=os(""),l=os([]);function o(){try{l.value=Ft()||[]}catch(i){console.error("Failed to load category data:",i),l.value=[]}}function p(i){return Qn(`/article/${i.replace(/\.md$/,"")}`)}function c(i){const d=t.value.replace(/\.md$/,""),h=i.replace(/\.md$/,""),y=d.replace(/-en$/,""),j=h.replace(/-en$/,"");return y===j}function u(){const i=t.value;if(!i){e.value=[];return}const d=i.split("/").filter(Boolean);if(d.length<2){e.value=[];return}const h=d[0],y=d[1];let j=null;if(Array.isArray(l.value)){s:for(const w of l.value)if(Array.isArray(w.items))for(const f of w.items){if(f?.name!==y)continue;let x=!1;if(Array.isArray(f.articles)&&(x=f.articles.some(k=>typeof k?.articleUrl=="string"&&k.articleUrl.includes(`/article/${h}/`))),!x&&Array.isArray(f.categories)&&(x=f.categories.some(k=>Array.isArray(k.articles)&&k.articles.some(_=>typeof _?.articleUrl=="string"&&_.articleUrl.includes(`/article/${h}/${y}/`)))),x){j=f;break s}}}if(!j){e.value=[];return}const T=[],g=[],m=(w,f)=>{const x=String(f).replace(/^\/+/,"").split("/"),k=x[0]==="article"?1:0,_=x[k],D=x[k+1],S=x.slice(k+2),I=`${_}/${D}/${S.join("/")}`;return{title:w,path:`${I}.md`}};Array.isArray(j.articles)&&j.articles.forEach(w=>{w?.articleUrl&&T.push(m(w.title,w.articleUrl))}),Array.isArray(j.categories)&&j.categories.forEach(w=>{const f=[];Array.isArray(w.articles)&&w.articles.forEach(x=>{x?.articleUrl&&f.push(m(x.title,x.articleUrl))}),f.length&&g.push({name:w.title||w.key,type:"directory",files:f})}),e.value=[{name:j.title||j.name||y,files:T,children:g}]}function r(){const i=a.params.path;t.value=Array.isArray(i)?i.join("/"):i||"",u()}return Hs(()=>a.params.path,r,{immediate:!0}),Hs(n,()=>{o(),r()}),wn(()=>{o(),u()}),(i,d)=>{const h=Le("router-link");return N(),B("div",t_,[(N(!0),B(bs,null,Ms(e.value,y=>(N(),B("div",{key:y.name,class:"tree-group"},[E("div",l_,X(y.name),1),E("ul",o_,[(N(!0),B(bs,null,Ms(y.children,j=>(N(),B(bs,{key:j.name},[j.files&&j.files.length?(N(),B("li",p_,[E("div",c_,[E("div",r_,X(j.name),1),j.files&&j.files.length?(N(),B("ul",i_,[(N(!0),B(bs,null,Ms(j.files,T=>(N(),B("li",{key:T.path,class:"article-item"},[ws(h,{to:p(T.path),class:Us(["tree-item tree-item--child",{"tree-item--active":c(T.path)}])},{default:vn(()=>[As(X(T.title),1)]),_:2},1032,["to","class"])]))),128))])):hs("",!0),j.children&&j.children.length?(N(),B("ul",u_,[(N(!0),B(bs,null,Ms(j.children,T=>(N(),B("li",{key:T.name,class:"article-item"},[E("div",h_,[E("div",d_,X(T.name),1),T.files&&T.files.length?(N(),B("ul",m_,[(N(!0),B(bs,null,Ms(T.files,g=>(N(),B("li",{key:g.path,class:"article-item"},[ws(h,{to:p(g.path),class:Us(["tree-item tree-item--child",{"tree-item--active":c(g.path)}])},{default:vn(()=>[As(X(g.title),1)]),_:2},1032,["to","class"])]))),128))])):hs("",!0)])]))),128))])):hs("",!0)])])):hs("",!0)],64))),128)),(N(!0),B(bs,null,Ms(y.files,j=>(N(),B("li",{key:j.path,class:"article-item"},[ws(h,{to:p(j.path),class:Us(["tree-item tree-item--child",{"tree-item--active":c(j.path)}])},{default:vn(()=>[As(X(j.title),1)]),_:2},1032,["to","class"])]))),128))])]))),128))])}}}),Cn=(s,n)=>{const a=s.__vccOpts||s;for(const[e,t]of n)a[e]=t;return a},Xi=Cn(j_,[["__scopeId","data-v-7ae1de2a"]]),f_={class:"app-header"},g_={class:"container px-0"},b_={class:"navbar navbar-expand-lg app-nav"},__={class:"container-fluid d-flex app-nav__inner"},y_={class:"app-nav__wordmark"},v_={class:"d-flex d-lg-none ms-auto app-nav__actions app-nav__actions--mobile"},w_=["aria-label"],C_=["aria-label"],k_={key:0,class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},x_={key:1,class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},P_=["aria-label"],E_=["aria-label"],S_={class:"navbar-nav mb-2 mb-lg-0 app-nav__links"},T_={key:0,class:"app-nav__link-divider","aria-hidden":"true"},A_={class:"d-none d-lg-flex ms-auto app-nav__actions"},R_=["aria-label"],L_=["aria-label"],M_={key:0,class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},D_={key:1,class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},O_=["aria-label"],I_={class:"offcanvas-panel"},N_={class:"offcanvas-section"},F_={class:"offcanvas-head"},$_={class:"offcanvas-card"},B_={class:"list-unstyled m-0"},q_={key:0,class:"offcanvas-section"},z_={class:"offcanvas-head"},V_=Gs({__name:"AppHeader",emits:["open-search"],setup(s,{emit:n}){const a=Ka(),e=Ie(),{t}=nn(),l=wo(),o=ss(()=>l.theme==="dark"?!0:l.theme==="light"?!1:typeof window<"u"&&window.matchMedia("(prefers-color-scheme: dark)").matches),p=os(!1),c=os(!1),u=n,r=k=>{k.currentTarget?.blur()},i=ss(()=>[{text:t("categories"),href:Qn("/category")},{text:t("resources"),href:Qn("/resource")},{text:t("about"),href:Qn("/about")}]),d=k=>k===a.path?!0:k!==Qn("/")&&a.path.startsWith(k),h=ss(()=>a.path.includes("/article/")),y=()=>{l.toggleTheme()},j=()=>{const k=a.path.match(/^\/(zh|en)/),_=k&&k[1]==="en"?"/zh":"/en",D=k?a.path.slice(k[0].length):a.path;e.push({path:`${_}${D}`,query:a.query})},T=()=>{window.innerWidth<992?w():p.value=!p.value},g=()=>{const k=document.documentElement,_=document.body;k&&(k.style.overflow="hidden",k.style.overscrollBehavior="contain"),_&&(_.style.overflow="hidden",_.style.overscrollBehavior="contain")},m=()=>{const k=document.documentElement,_=document.body;k&&(k.style.overflow="",k.style.overscrollBehavior=""),_&&(_.style.overflow="",_.style.overscrollBehavior="")},w=()=>{g(),c.value=!0},f=()=>{c.value=!1,m()},x=k=>{const _=k.target;_&&_.closest("a")&&f()};return wn(()=>{l.initTheme(),l.initLocale()}),qn(()=>{m()}),(k,_)=>{const D=Le("RouterLink");return N(),B(bs,null,[E("header",f_,[E("div",g_,[E("nav",b_,[E("div",__,[ws(D,{class:"navbar-brand app-nav__brand",to:cs(Qn)("/"),onClick:_[0]||(_[0]=S=>p.value=!1)},{default:vn(()=>[E("span",y_,[As(X(cs(un).author),1),_[4]||(_[4]=E("span",{class:"app-nav__apos"},"’",-1)),_[5]||(_[5]=As("s blog",-1))])]),_:1},8,["to"]),E("div",v_,[E("button",{class:"icon-btn",onClick:_[1]||(_[1]=S=>u("open-search")),onFocus:r,"aria-label":cs(t)("search")},[..._[6]||(_[6]=[E("svg",{class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[E("circle",{cx:"11",cy:"11",r:"8"}),E("line",{x1:"21",y1:"21",x2:"16.65",y2:"16.65"})],-1)])],40,w_),E("button",{class:"icon-btn",onClick:y,onFocus:r,"aria-label":cs(t)("theme")},[o.value?(N(),B("svg",k_,[..._[7]||(_[7]=[He('<circle cx="12" cy="12" r="5" data-v-b2e68a7c></circle><line x1="12" y1="1" x2="12" y2="3" data-v-b2e68a7c></line><line x1="12" y1="21" x2="12" y2="23" data-v-b2e68a7c></line><line x1="4.22" y1="4.22" x2="5.64" y2="5.64" data-v-b2e68a7c></line><line x1="18.36" y1="18.36" x2="19.78" y2="19.78" data-v-b2e68a7c></line><line x1="1" y1="12" x2="3" y2="12" data-v-b2e68a7c></line><line x1="21" y1="12" x2="23" y2="12" data-v-b2e68a7c></line><line x1="4.22" y1="19.78" x2="5.64" y2="18.36" data-v-b2e68a7c></line><line x1="18.36" y1="5.64" x2="19.78" y2="4.22" data-v-b2e68a7c></line>',9)])])):(N(),B("svg",x_,[..._[8]||(_[8]=[E("path",{d:"M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"},null,-1)])]))],40,C_),E("button",{class:"icon-btn",onClick:j,onFocus:r,"aria-label":cs(t)("language")},[..._[9]||(_[9]=[He('<svg class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.75" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" data-v-b2e68a7c><path d="m5 8 6 6" data-v-b2e68a7c></path><path d="m4 14 6-6 2-3" data-v-b2e68a7c></path><path d="M2 5h12" data-v-b2e68a7c></path><path d="M7 2h1" data-v-b2e68a7c></path><path d="m22 22-5-10-5 10" data-v-b2e68a7c></path><path d="M14 18h6" data-v-b2e68a7c></path></svg>',1)])],40,P_),E("button",{class:"icon-btn",onClick:T,onFocus:r,"aria-label":cs(t)("menu")},[..._[10]||(_[10]=[E("svg",{class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[E("line",{x1:"3",y1:"12",x2:"21",y2:"12"}),E("line",{x1:"3",y1:"6",x2:"21",y2:"6"}),E("line",{x1:"3",y1:"18",x2:"21",y2:"18"})],-1)])],40,E_)]),E("div",{class:Us(["navbar-collapse collapse",{show:p.value}])},[E("ul",S_,[i.value.length?(N(),B("li",T_)):hs("",!0),(N(!0),B(bs,null,Ms(i.value,S=>(N(),B("li",{class:"nav-item",key:S.text},[ws(D,{class:Us(["nav-link app-nav__link",{active:d(S.href)}]),to:S.href,onClick:_[2]||(_[2]=I=>p.value=!1)},{default:vn(()=>[As(X(S.text),1)]),_:2},1032,["to","class"])]))),128))])],2),E("div",A_,[E("button",{class:"icon-btn",onClick:_[3]||(_[3]=S=>u("open-search")),onFocus:r,"aria-label":cs(t)("search")},[..._[11]||(_[11]=[E("svg",{class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[E("circle",{cx:"11",cy:"11",r:"8"}),E("line",{x1:"21",y1:"21",x2:"16.65",y2:"16.65"})],-1)])],40,R_),E("button",{class:"icon-btn",onClick:y,onFocus:r,"aria-label":cs(t)("theme")},[o.value?(N(),B("svg",M_,[..._[12]||(_[12]=[He('<circle cx="12" cy="12" r="5" data-v-b2e68a7c></circle><line x1="12" y1="1" x2="12" y2="3" data-v-b2e68a7c></line><line x1="12" y1="21" x2="12" y2="23" data-v-b2e68a7c></line><line x1="4.22" y1="4.22" x2="5.64" y2="5.64" data-v-b2e68a7c></line><line x1="18.36" y1="18.36" x2="19.78" y2="19.78" data-v-b2e68a7c></line><line x1="1" y1="12" x2="3" y2="12" data-v-b2e68a7c></line><line x1="21" y1="12" x2="23" y2="12" data-v-b2e68a7c></line><line x1="4.22" y1="19.78" x2="5.64" y2="18.36" data-v-b2e68a7c></line><line x1="18.36" y1="5.64" x2="19.78" y2="4.22" data-v-b2e68a7c></line>',9)])])):(N(),B("svg",D_,[..._[13]||(_[13]=[E("path",{d:"M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"},null,-1)])]))],40,L_),E("button",{class:"icon-btn",onClick:j,onFocus:r,"aria-label":cs(t)("language")},[..._[14]||(_[14]=[He('<svg class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.75" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" data-v-b2e68a7c><path d="m5 8 6 6" data-v-b2e68a7c></path><path d="m4 14 6-6 2-3" data-v-b2e68a7c></path><path d="M2 5h12" data-v-b2e68a7c></path><path d="M7 2h1" data-v-b2e68a7c></path><path d="m22 22-5-10-5 10" data-v-b2e68a7c></path><path d="M14 18h6" data-v-b2e68a7c></path></svg>',1)])],40,O_)])])])])]),c.value?(N(),B("div",{key:0,class:"mobile-offcanvas d-lg-none",onClick:_n(f,["self"])},[E("div",I_,[E("div",N_,[E("div",F_,X(cs(t)("menu")),1),E("div",$_,[E("ul",B_,[(N(!0),B(bs,null,Ms(i.value,S=>(N(),B("li",{key:S.text,class:"my-1"},[ws(D,{to:S.href,class:Us(["offcanvas-link d-flex align-items-center",{active:d(S.href)}]),onClick:f},{default:vn(()=>[E("span",null,X(S.text),1),_[15]||(_[15]=E("i",{class:"fas fa-chevron-right offcanvas-link__chevron"},null,-1))]),_:2},1032,["to","class"])]))),128))])])]),h.value?(N(),B("div",q_,[E("div",z_,X(cs(t)("tableOfContents")),1),E("div",{class:"offcanvas-tree offcanvas-card",onClick:x},[ws(Xi)])])):hs("",!0)]),E("div",{class:"offcanvas-backdrop",onClick:f})])):hs("",!0)],64)}}}),U_=Cn(V_,[["__scopeId","data-v-b2e68a7c"]]),H_={class:"site-footer"},G_={class:"container"},W_={class:"site-footer__inner"},K_={class:"footer-copy"},X_={class:"footer-designed"},Y_={class:"footer-designed__text"},Q_={class:"footer-designed__name"},J_={class:"footer-designed__icons"},Z_=["href"],sy=["href"],ny=Gs({__name:"AppFooter",setup(s){const{t:n}=nn(),a=new Date().getFullYear(),e=un.startYear&&un.startYear<a?`${un.startYear} - ${a}`:`${a}`;return(t,l)=>(N(),B("footer",H_,[E("div",G_,[E("div",W_,[E("p",K_,"© "+X(cs(e))+" "+X(cs(un).author),1),E("div",X_,[E("p",Y_,[As(X(cs(n)("designedByPrefix")),1),E("strong",Q_,X(cs(un).author),1),As(X(cs(n)("designedBySuffix")),1)]),E("div",J_,[E("a",{href:cs(un).github,target:"_blank",rel:"noopener noreferrer",class:"footer-link","aria-label":"GitHub"},[...l[0]||(l[0]=[E("svg",{class:"footer-icon",width:"16",height:"16",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[E("path",{d:"M9 19c-5 1.5-5-2.5-7-3m14 6v-3.87a3.37 3.37 0 0 0-.94-2.61c3.14-.35 6.44-1.54 6.44-7A5.44 5.44 0 0 0 20 4.77 5.07 5.07 0 0 0 19.91 1S18.73.65 16 2.48a13.38 13.38 0 0 0-7 0C6.27.65 5.09 1 5.09 1A5.07 5.07 0 0 0 5 4.77a5.44 5.44 0 0 0-1.5 3.78c0 5.42 3.3 6.61 6.44 7A3.37 3.37 0 0 0 9 18.13V22"})],-1)])],8,Z_),E("a",{href:`mailto:${cs(un).email}`,class:"footer-link","aria-label":"Email"},[...l[1]||(l[1]=[E("svg",{class:"footer-icon",width:"16",height:"16",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[E("path",{d:"M4 4h16c1.1 0 2 .9 2 2v12c0 1.1-.9 2-2 2H4c-1.1 0-2-.9-2-2V6c0-1.1.9-2 2-2z"}),E("polyline",{points:"22,6 12,13 2,6"})],-1)])],8,sy)])])])])]))}}),ay=Cn(ny,[["__scopeId","data-v-f017e1f3"]]),ey=["aria-label"],fc="btt",ty=Gs({__name:"BackToTop",setup(s){const{t:n}=nn(),a=os(!1),e=os(null),t=os(!1),l=os(!1),o=os(0),p=os(0),c=os((typeof window<"u"?window.innerHeight:1024)-100),u=os(!1),r=os(null),i=os(null);function d(){return{gap:48,minTop:68,maxTop:window.innerHeight-40-20}}function h(k){const{minTop:_,maxTop:D}=d();return Math.max(_,Math.min(D,k))}function y(k){e.value=k,!a.value&&(a.value=!0,requestAnimationFrame(()=>{window.dispatchEvent(new CustomEvent("floating-buttons-base-top",{detail:{baseTop:e.value,source:fc}})),a.value=!1}))}function j(k){const _=k.detail,D=_?.baseTop;_?.source!==fc&&typeof D=="number"&&(c.value=h(D))}function T(){if(r.value)return;const k=document.createElement("div");k.style.cssText="position:absolute;left:0;top:0;width:1px;height:181px;pointer-events:none;visibility:hidden",document.body.appendChild(k),i.value=k,r.value=new IntersectionObserver(([_])=>{t.value=!_.isIntersecting,t.value&&y(c.value)},{threshold:0}),r.value.observe(k)}function g(){const k=typeof window.matchMedia=="function"&&window.matchMedia("(prefers-reduced-motion: reduce)").matches;window.scrollTo({top:0,behavior:k?"auto":"smooth"})}function m(){l.value||g()}function w(k){k.preventDefault(),l.value=!0,u.value=!1,o.value=k.touches[0].clientY,p.value=c.value}function f(k){u.value=!0;const D=k.touches[0].clientY-o.value;c.value=h(p.value+D),y(c.value),k.preventDefault()}function x(k){k.preventDefault(),l.value=!1,u.value||g()}return wn(()=>{T(),window.addEventListener("floating-buttons-base-top",j),c.value=window.innerHeight-100,y(c.value)}),qn(()=>{r.value&&(r.value.disconnect(),r.value=null),i.value&&(i.value.remove(),i.value=null),window.removeEventListener("floating-buttons-base-top",j)}),(k,_)=>Ua((N(),B("button",{class:"back-to-top d-flex align-items-center justify-content-center",onClick:m,"aria-label":cs(n)("backToTop"),onTouchstart:_n(w,["prevent","stop"]),onTouchmove:_n(f,["prevent","stop"]),onTouchend:_n(x,["prevent","stop"]),style:ba({top:c.value+"px"})},[..._[0]||(_[0]=[E("i",{class:"fas fa-arrow-up"},null,-1)])],44,ey)),[[qr,t.value]])}}),ly=Cn(ty,[["__scopeId","data-v-7a6aaf37"]]),oy={class:"page-wrap"},py={class:"main-content"},cy=Gs({__name:"App",setup(s){const n=os(!1),a=gh(()=>ys(()=>import("./SearchModal-CPATOImV.js"),__vite__mapDeps([0,1]))),e=Ka(),t=wo(),{locale:l}=nn(),o=ss(()=>{const i=e.path.replace(/^\/(zh|en)(?=\/|$)/,"");return i==="/"||i===""?"":i}),p=ss(()=>l.value==="en-US"?"/en":"/zh");Tt({htmlAttrs:{lang:()=>l.value==="en-US"?"en":"zh-CN"},link:ss(()=>{const i=o.value.includes("/article/"),d=i?o.value.replace(/-en$/,""):o.value,h=i&&!o.value.endsWith("-en")?`${o.value}-en`:o.value;return[{rel:"canonical",href:`${un.url}${p.value}${o.value}`},{rel:"alternate",hreflang:"zh-CN",href:`${un.url}/zh${d}`},{rel:"alternate",hreflang:"en",href:`${un.url}/en${h}`},{rel:"alternate",hreflang:"x-default",href:`${un.url}/zh${d}`}]}),meta:ss(()=>[{property:"og:url",content:`${un.url}${p.value}${o.value}`}])});const c=i=>{const d=i.match(/^\/(zh|en)(\/|$)/);d&&t.setLocale($i(d[1]))},u=i=>i instanceof HTMLElement?i.tagName==="INPUT"||i.tagName==="TEXTAREA"||i.isContentEditable:!1,r=i=>{i.key==="/"&&!u(i.target)&&(i.preventDefault(),n.value=!0)};return Hs(()=>e.fullPath,c),wn(()=>{c(e.fullPath),window.addEventListener("keydown",r)}),qn(()=>{window.removeEventListener("keydown",r)}),(i,d)=>{const h=Le("router-view");return N(),B("div",oy,[ws(U_,{onOpenSearch:d[0]||(d[0]=y=>n.value=!0)}),E("main",py,[ws(h,null,{default:vn(({Component:y})=>[ws(Cd,{name:"view",mode:"out-in"},{default:vn(()=>[(N(),Sa(Eh(y)))]),_:2},1024)]),_:1})]),ws(ay),ws(ly),n.value?(N(),Sa(cs(a),{key:0,onClose:d[1]||(d[1]=y=>n.value=!1)})):hs("",!0)])}}}),ry=Cn(cy,[["__scopeId","data-v-73d8966b"]]);function Yi(s){return Math.max(1,Math.round(Number(s)/300))}const iy={class:"post-item__index num","aria-hidden":"true"},uy={class:"post-item__meta"},hy={class:"post-item__cat"},dy={class:"post-item__date"},my={class:"post-item__words"},jy={class:"post-item__reading"},fy={class:"post-item__title"},gy={class:"post-item__preview"},by={class:"post-item__tags"},_y=["onClick"],yy=["aria-label"],vy={class:"pagination__side"},wy=["aria-label"],Cy={class:"pagination__pages"},ky={key:1,class:"page-ellipsis"},xy=["onClick"],Py={key:2,class:"page-ellipsis"},Ey={class:"pagination__side pagination__side--right"},Sy=["aria-label"],Ty=Gs({__name:"PostList",props:{docs:{},perPage:{default:6}},setup(s){const n=s,{t:a,locale:e}=nn(),t=Ka(),l=Ie(),o=os(1),p=os(5),c=os([]),u=os([]),r=ss(()=>a("pagination")),i=ss(()=>Math.max(1,Math.ceil(n.docs.length/n.perPage))),d=ss(()=>{const F=(o.value-1)*n.perPage,G=F+n.perPage;return n.docs.slice(F,G)}),h=ss(()=>{const F=[],G=i.value,is=o.value,ps=p.value;if(G<=ps){for(let us=1;us<=G;us++)F.push(us);return F}F.push(is);let Y=1;for(;F.length<ps&&(is-Y>=1&&!F.includes(is-Y)&&F.push(is-Y),!(F.length>=ps));)is+Y<=G&&!F.includes(is+Y)&&F.push(is+Y),Y++;return F.includes(1)||(F.pop(),F.push(1)),F.length<ps&&!F.includes(G)?F.push(G):F.length>=ps&&!F.includes(G)&&(F.pop(),F.push(G)),F.sort((us,Ns)=>us-Ns)}),y=ss(()=>h.value.filter(F=>F!==1&&F!==i.value)),j=ss(()=>h.value.includes(1)),T=ss(()=>h.value.includes(i.value)),g=ss(()=>i.value>p.value&&h.value[1]>2),m=ss(()=>i.value>p.value&&h.value[h.value.length-2]<i.value-1),w=ss(()=>{const F={};try{(u.value||[]).forEach(G=>{(G.items||[]).forEach(is=>{(is.categories||[]).forEach(ps=>{ps&&ps.key&&ps.title&&(F[ps.key]=ps.title)})})})}catch(G){console.error(G)}return F});function f(){try{c.value=Zb()||[],u.value=Ft()||[]}catch(F){console.error("Failed to load data:",F),c.value=[],u.value=[]}}function x(F){let G="";const is=Array.isArray(c.value)?c.value.find(ps=>ps.title===F.title):null;if(is&&is.relativePath)G=`notes/${is.relativePath}.md`;else{const ps=e.value==="en-US",Y=F.category?.[1]||"notes",us=F.title.toLowerCase().replace(/[^a-z0-9]/g,"-");G=`${Y}/${us}.md`,ps&&(G=G.replace(".md","-en.md"))}return Qn(`/article/${G.replace(/\.md$/,"")}`)}function k(F){if(F>=1&&F<=i.value){o.value=F;const G={...t.query,page:String(F)};l.push({path:t.path,query:G}).catch(()=>{}),Js(()=>window.scrollTo({top:0,behavior:"smooth"}))}}function _(){if(o.value>1){const F=o.value-1;o.value=F;const G={...t.query,page:String(F)};l.push({path:t.path,query:G}).catch(()=>{}),Js(()=>window.scrollTo({top:0,behavior:"smooth"}))}}function D(){if(o.value<i.value){const F=o.value+1;o.value=F;const G={...t.query,page:String(F)};l.push({path:t.path,query:G}).catch(()=>{}),Js(()=>window.scrollTo({top:0,behavior:"smooth"}))}}function S(F){if(!F)return;const G={...t.query,tag:F,page:"1"},is=e.value==="en-US"?"/en":"/zh";l.push({path:`${is}/`,query:G}).catch(()=>{}),Js(()=>window.scrollTo({top:0,behavior:"smooth"}))}function I(){p.value=window.innerWidth<480?3:5}function J(F){if(!Array.isArray(F)||F.length===0)return"";const G=F[0]||"",is=F[1];if(!is)return G;const ps=w.value[is]||is;return`${G} / ${ps}`}function U(F){return a("postReadingTime",{minutes:Yi(F)})}return f(),Hs(()=>n.docs,()=>{o.value=1}),Hs(()=>t.query.page,F=>{const G=parseInt(String(F)),is=Number.isFinite(G)&&G>=1?Math.min(G,i.value):1;is!==o.value&&(o.value=is,Js(()=>window.scrollTo({top:0,behavior:"smooth"})))}),Hs(e,()=>f()),wn(()=>{const F=parseInt(String(t.query.page));o.value=Number.isFinite(F)&&F>=1?Math.min(F,i.value):1,I(),window.addEventListener("resize",I)}),qn(()=>{window.removeEventListener("resize",I)}),(F,G)=>{const is=Le("router-link"),ps=lo("reveal");return N(),B("div",null,[(N(!0),B(bs,null,Ms(d.value,(Y,us)=>Ua((N(),B("article",{key:Y.id,class:"post-item",style:ba({"--reveal-delay":Math.min(us,5)*40+"ms"})},[E("span",iy,X(String((o.value-1)*n.perPage+us+1).padStart(2,"0")),1),ws(is,{to:x(Y),class:"post-item__main"},{default:vn(()=>[E("div",uy,[E("span",hy,X(J(Y.category)),1),G[5]||(G[5]=E("span",{class:"divider-v"},null,-1)),E("span",dy,[G[2]||(G[2]=E("i",{class:"fas fa-calendar-alt"},null,-1)),As(X(Y.date),1)]),G[6]||(G[6]=E("span",{class:"divider-v"},null,-1)),E("span",my,[G[3]||(G[3]=E("i",{class:"fas fa-file-lines"},null,-1)),As(X(Y.wordCount)+" "+X(cs(a)("wordUnit")),1)]),G[7]||(G[7]=E("span",{class:"divider-v"},null,-1)),E("span",jy,[G[4]||(G[4]=E("i",{class:"fas fa-clock"},null,-1)),As(X(U(Y.wordCount)),1)])]),E("h3",fy,X(Y.title),1),E("p",gy,X(Y.preview),1)]),_:2},1032,["to"]),E("div",by,[(N(!0),B(bs,null,Ms(Y.tags,Ns=>(N(),B("span",{key:Ns,class:"post-item__tag",onClick:Xs=>S(Ns)},X(Ns),9,_y))),128))]),ws(is,{to:x(Y),class:"post-item__arrow","aria-label":Y.title},{default:vn(()=>[...G[8]||(G[8]=[E("i",{class:"fas fa-arrow-right"},null,-1)])]),_:1},8,["to","aria-label"])],4)),[[ps]])),128)),i.value>1?(N(),B("nav",{key:0,class:"pagination","aria-label":r.value},[E("div",vy,[o.value>1?(N(),B("button",{key:0,class:"page-btn page-btn--nav",onClick:_,"aria-label":cs(a)("prevPage")},[...G[9]||(G[9]=[E("i",{class:"fas fa-chevron-left"},null,-1)])],8,wy)):hs("",!0)]),E("div",Cy,[j.value?(N(),B("button",{key:0,class:Us(["page-btn",{"page-btn--active":o.value===1}]),onClick:G[0]||(G[0]=Y=>k(1))},"1",2)):hs("",!0),g.value?(N(),B("span",ky,"...")):hs("",!0),(N(!0),B(bs,null,Ms(y.value,Y=>(N(),B("button",{key:Y,class:Us(["page-btn",{"page-btn--active":o.value===Y}]),onClick:us=>k(Y)},X(Y),11,xy))),128)),m.value?(N(),B("span",Py,"...")):hs("",!0),T.value&&i.value>1?(N(),B("button",{key:3,class:Us(["page-btn",{"page-btn--active":o.value===i.value}]),onClick:G[1]||(G[1]=Y=>k(i.value))},X(i.value),3)):hs("",!0)]),E("div",Ey,[o.value<i.value?(N(),B("button",{key:0,class:"page-btn page-btn--nav",onClick:D,"aria-label":cs(a)("nextPage")},[...G[10]||(G[10]=[E("i",{class:"fas fa-chevron-right"},null,-1)])],8,Sy)):hs("",!0)])],8,yy)):hs("",!0)])}}}),Ay=Cn(Ty,[["__scopeId","data-v-1043b469"]]),Ry={class:"page-section home-section"},Ly={class:"hero"},My={class:"hero__main"},Dy={class:"hero__greeting"},Oy={class:"hero__greeting-mark"},Iy={class:"hero__name"},Ny={class:"hero__bio"},Fy={class:"hero__stats"},$y={class:"hero__stat"},By={class:"hero__stat-num num"},qy={class:"hero__stat-label"},zy={class:"hero__stat"},Vy={class:"hero__stat-num num"},Uy={class:"hero__stat-label"},Hy={class:"hero__stat"},Gy={class:"hero__stat-num num"},Wy={class:"hero__stat-label"},Ky={key:0,class:"tags-row"},Xy=["onClick"],Yy={class:"posts-header"},Qy={key:0,class:"posts-header__tag-title"},Jy={class:"posts-header__actions"},Zy=["aria-label"],sv=["aria-label"],nv=Gs({name:"HomeView",__name:"Home",setup(s){Tt({title:"zorrooz’s blog - Home"});const{t:n,locale:a}=nn(),e=Ka(),t=Ie(),l=os([]),o=os([]),p=un.author,c=ss(()=>e.query.tag||""),u=ss(()=>typeof e.query.from=="string"?e.query.from:""),r=ss(()=>{const f=c.value;return f?l.value.filter(x=>Array.isArray(x.tags)&&x.tags.includes(f)):l.value}),i=ss(()=>{const f=new Map;l.value.forEach(S=>(S.tags||[]).forEach(I=>f.set(I,(f.get(I)||0)+1)));const x=Array.from(f.entries()).map(([S,I])=>({name:S,count:I})).sort((S,I)=>S.name.localeCompare(I.name)),k=x.length,_=Math.ceil(k/3),D=new Map([...x].sort((S,I)=>I.count-S.count).map((S,I)=>[S.name,I]));return x.map(S=>{const I=D.get(S.name)??0;return{...S,level:I<_?"lg":I<_*2?"md":"sm"}})}),d=ss(()=>Array.isArray(l.value)?l.value.length:0),h=ss(()=>Array.isArray(o.value)?o.value.length:0),y=ss(()=>Array.isArray(l.value)?l.value.reduce((f,x)=>{const k=typeof x?.wordCount=="number"?x.wordCount:0;return f+(Number.isFinite(k)?k:0)},0):0),j=ss(()=>{const f=y.value;return f>=1e6?(f/1e6).toFixed(f%1e6?1:0)+"M":f>=1e3?(f/1e3).toFixed(f%1e3?1:0)+"K":String(f)});function T(){try{l.value=Jb()||[],o.value=s_()||[]}catch(f){console.error("Failed to load post data:",f),l.value=[],o.value=[]}}function g(){const f={...e.query};delete f.tag,delete f.from,f.page="1",t.push({path:e.path,query:f}).catch(()=>{}),Js(()=>{window.scrollTo({top:0,behavior:"smooth"}),T()})}function m(){u.value&&t.push(u.value).catch(()=>{})}function w(f){if(!f)return;const x={...e.query,tag:f,page:"1"},k=a.value==="en-US"?"/en":"/zh";t.push({path:`${k}/`,query:x}).catch(()=>{}),Js(()=>window.scrollTo({top:0,behavior:"smooth"}))}return T(),Hs(a,(f,x)=>{f!==x&&(c.value?g():T())}),wn(()=>{T()}),(f,x)=>{const k=lo("reveal");return N(),B("div",Ry,[E("header",Ly,[E("div",My,[E("span",Dy,[E("span",Oy,X(cs(n)("greetingPrefix")),1),As(X(cs(n)("greeting")),1)]),E("h1",Iy,X(cs(p)),1),E("p",Ny,X(cs(n)("developer")),1)]),E("div",Fy,[E("div",$y,[E("span",By,X(d.value),1),E("span",qy,X(cs(n)("articles")),1)]),E("div",zy,[E("span",Vy,X(h.value),1),E("span",Uy,X(cs(n)("tags")),1)]),E("div",Hy,[E("span",Gy,X(j.value),1),E("span",Wy,X(cs(n)("words")),1)])])]),c.value?hs("",!0):Ua((N(),B("div",Ky,[(N(!0),B(bs,null,Ms(i.value,_=>(N(),B("span",{key:_.name,class:Us(["tag",`tag--${_.level}`]),onClick:D=>w(_.name)},X(_.name),11,Xy))),128))])),[[k]]),Ua((N(),B("div",Yy,[E("h2",{class:Us(["posts-header__title",{"posts-header__title--tag":c.value}])},[c.value?(N(),B("span",Qy,"# "+X(c.value),1)):(N(),B(bs,{key:1},[As(X(cs(n)("recentPosts")),1)],64))],2),E("div",Jy,[u.value?(N(),B("button",{key:0,class:"chip-close",onClick:m,"aria-label":cs(n)("backToArticle")},[x[0]||(x[0]=E("i",{class:"fas fa-arrow-left"},null,-1)),As(X(cs(n)("backToArticle")),1)],8,Zy)):hs("",!0),c.value?(N(),B("button",{key:1,class:"chip-close",onClick:g,"aria-label":cs(n)("close")},[x[1]||(x[1]=E("i",{class:"fas fa-times"},null,-1)),As(X(cs(n)("clearFilter")),1)],8,sv)):hs("",!0)])])),[[k]]),ws(Ay,{docs:r.value,perPage:5},null,8,["docs"])])}}}),av=Cn(nv,[["__scopeId","data-v-f7cc79bf"]]),ev={class:"page-section category-view"},tv={class:"category-head"},lv={class:"article-title"},ov={class:"cat-section__header"},pv={class:"cat-section__title"},cv={class:"cat-section__count"},rv={class:"cat-grid"},iv={class:"cat-card__head"},uv={class:"cat-card__name"},hv={class:"cat-card__ext-links"},dv=["href"],mv=["href"],jv={class:"cat-card__desc"},fv={key:0,class:"cat-card__stats"},gv={key:0,class:"cat-stat"},bv={key:1,class:"cat-stat"},_v={key:2,class:"cat-stat"},yv={key:1,class:"cat-card__stats"},vv={key:0,class:"cat-stat"},wv={key:1,class:"cat-stat"},Cv={key:2,class:"cat-stat"},kv={key:2,class:"cat-card__tags"},xv={class:"cat-card__links"},Pv=["onClick"],Ev=Gs({name:"CategoryView",__name:"Category",setup(s){Tt({title:"zorrooz’s blog - Categories"});const{t:n,locale:a}=nn(),e=Ie(),t=os([]),l=ss(()=>n("categories")),o=ss(()=>n("notes")),p=ss(()=>n("projects")),c=ss(()=>n("topics")),u=ss(()=>n("seeMore"));function r(){try{const m=Ft();t.value=Array.isArray(m)?m:m?.categoryList||[]}catch(m){console.error("Failed to load category data:",m),t.value=[]}}function i(m){return m===o.value?"fa-book-open":m===p.value?"fa-folder-open":m===c.value?"fa-flask":"fa-folder"}function d(m){const w=m?.items||[];if(m.title===o.value){const f=w.reduce((x,k)=>x+(k.stats?.postsCount||0),0);return n("countPosts",{count:f})}return m.title===p.value?n("countProjects",{count:w.length}):n("countTopics",{count:w.length})}function h(m){return m&&m.stats&&m.stats.latestDate||""}function y(m){return typeof m!="string"||!m.trim()?"":/^https?:\/\//i.test(m)?m:"https://"+m.replace(/^\/+/,"")}function j(m){return typeof m!="string"||!m.trim()?"":/^https?:\/\//i.test(m)?m:/^10\.\d{4,9}\//.test(m)?"https://doi.org/"+m:"https://"+m.replace(/^\/+/,"")}function T(m){return!!(m?.url||m?.github||m?.doi)}function g(m){const w=m?.root;if(w){e.push(Qn(w)).catch(_=>{_.name!=="NavigationDuplicated"&&!_.toString().includes("Navigation cancelled")&&console.error("Navigation error:",_)});return}const f=y(m?.url);if(f){window.open(f,"_blank","noopener,noreferrer");return}const x=y(m?.github);if(x){window.open(x,"_blank","noopener,noreferrer");return}const k=j(m?.doi);k&&window.open(k,"_blank","noopener,noreferrer")}return r(),Hs(a,()=>{r()}),(m,w)=>{const f=lo("reveal");return N(),B("div",ev,[E("div",tv,[E("h1",lv,X(l.value),1)]),(N(!0),B(bs,null,Ms(t.value,(x,k)=>(N(),B("div",{key:k,class:"cat-section"},[E("div",ov,[E("h2",pv,[E("i",{class:Us([["fas",i(x.title)],"cat-section__icon"]),"aria-hidden":"true"},null,2),As(" "+X(x.title),1)]),E("span",cv,X(d(x)),1)]),E("div",rv,[(N(!0),B(bs,null,Ms(x.items,(_,D)=>Ua((N(),B("div",{key:D,class:"cat-card",style:ba({"--reveal-delay":Math.min(Number(D),5)*40+"ms"})},[E("div",iv,[E("div",uv,X(_.title),1),E("div",hv,[_.github?(N(),B("a",{key:0,href:y(_.github),target:"_blank",rel:"noopener noreferrer",class:"cat-card__ext-link","aria-label":"GitHub"},[...w[0]||(w[0]=[E("i",{class:"fab fa-github"},null,-1)])],8,dv)):hs("",!0),_.doi?(N(),B("a",{key:1,href:j(_.doi),target:"_blank",rel:"noopener noreferrer",class:"cat-card__ext-link","aria-label":"DOI"},[...w[1]||(w[1]=[E("i",{class:"fas fa-link"},null,-1)])],8,mv)):hs("",!0)])]),E("p",jv,X(_.desc),1),x.title===o.value?(N(),B("div",fv,[_.stats?.postsCount?(N(),B("span",gv,[w[2]||(w[2]=E("i",{class:"fas fa-file-lines"},null,-1)),As(X(cs(n)("countPosts",{count:_.stats.postsCount})),1)])):hs("",!0),_.stats?.totalWords?(N(),B("span",bv,[w[3]||(w[3]=E("i",{class:"fas fa-font"},null,-1)),As(X(cs(n)("countWords",{count:_.stats.totalWords})),1)])):hs("",!0),h(_)?(N(),B("span",_v,[w[4]||(w[4]=E("i",{class:"fas fa-clock"},null,-1)),As(X(h(_)),1)])):hs("",!0)])):(N(),B("div",yv,[_.language?(N(),B("span",vv,[w[5]||(w[5]=E("i",{class:"fas fa-code"},null,-1)),As(X(_.language),1)])):hs("",!0),_.year?(N(),B("span",wv,[w[6]||(w[6]=E("i",{class:"fas fa-calendar"},null,-1)),As(X(_.year),1)])):hs("",!0),_.license?(N(),B("span",Cv,[w[7]||(w[7]=E("i",{class:"fas fa-scale-balanced"},null,-1)),As(X(_.license),1)])):hs("",!0)])),Array.isArray(_.tags)&&_.tags.length?(N(),B("div",kv,[(N(!0),B(bs,null,Ms(_.tags,(S,I)=>(N(),B("span",{key:I,class:"cat-card__tag"},X(S),1))),128))])):hs("",!0),E("div",xv,[T(_)||_.root?(N(),B("a",{key:0,class:"cat-card__link",onClick:_n(S=>g(_),["prevent"])},[As(X(u.value),1),w[8]||(w[8]=E("i",{class:"fas fa-arrow-right"},null,-1))],8,Pv)):hs("",!0)])],4)),[[f]])),128))])]))),128))])}}}),Sv=Cn(Ev,[["__scopeId","data-v-fd52f9cf"]]),Tv={class:"page-section resource-view"},Av={class:"resource-head"},Rv={class:"article-title"},Lv={class:"resource-subtitle"},Mv={class:"res-layout"},Dv={class:"res-sidebar"},Ov={class:"res-group__label"},Iv={class:"res-group__count num"},Nv={class:"res-group__items"},Fv=["onClick"],$v={class:"res-main"},Bv={class:"res-groups"},qv={class:"res-group-block__title"},zv={class:"res-grid"},Vv=["href"],Uv={class:"res-card__head"},Hv={class:"res-card__name"},Gv={class:"res-card__ext-links"},Wv={key:0,class:"res-card__ext-link","aria-label":"DOI"},Kv={class:"res-card__desc"},Xv={class:"res-card__footer"},Yv={class:"res-card__url"},Qv=Gs({name:"ResourceView",__name:"Resource",setup(s){const{t:n,locale:a}=nn(),e=os([]),t=os(null),l=ss(()=>n("resources")),o=ss(()=>n("resourceSubtitle")),p=ss(()=>t.value?Array.isArray(t.value.children)&&t.value.children.length?t.value.children.filter(h=>Array.isArray(h.items)&&h.items.length):Array.isArray(t.value.items)?[{title:t.value.title,items:t.value.items}]:[]:[]);function c(h){t.value=h}function u(h){return t.value===h}function r(h){return h?h.replace(/^https?:\/\//,"").replace(/\/$/,""):""}function i(h){return!!h&&(h.includes("doi.org")||/^10\.\d{4,9}\//.test(h))}function d(){try{e.value=a_()||[],t.value=e.value[0]?.children?.[0]||null}catch(h){console.error("Failed to load resources data:",h),e.value=[],t.value=null}}return d(),Hs(a,()=>{d()}),(h,y)=>(N(),B("div",Tv,[E("header",Av,[E("h1",Rv,X(l.value),1),E("p",Lv,[y[0]||(y[0]=E("i",{class:"fas fa-circle-info resource-head__icon"},null,-1)),As(X(o.value),1)])]),E("div",Mv,[E("aside",Dv,[(N(!0),B(bs,null,Ms(e.value,j=>(N(),B("div",{key:j.title,class:"res-group"},[E("div",Ov,[E("span",null,X(j.title),1),E("span",Iv,X(j.children?.length||0),1)]),E("div",Nv,[(N(!0),B(bs,null,Ms(j.children,T=>(N(),B("button",{key:T.title,class:Us(["res-item",{"res-item--active":u(T)}]),onClick:g=>c(T)},X(T.title),11,Fv))),128))])]))),128))]),y[3]||(y[3]=E("div",{class:"res-divider","aria-hidden":"true"},null,-1)),E("main",$v,[E("div",Bv,[(N(!0),B(bs,null,Ms(p.value,j=>(N(),B("section",{key:j.title,class:"res-group-block"},[E("h3",qv,X(j.title),1),E("div",zv,[(N(!0),B(bs,null,Ms(j.items,T=>(N(),B("a",{key:T.name,href:T.url,target:"_blank",rel:"noopener noreferrer",class:"res-card"},[E("div",Uv,[E("span",Hv,X(T.name),1),E("span",Gv,[i(T.url)?(N(),B("span",Wv,[...y[1]||(y[1]=[E("i",{class:"fas fa-link"},null,-1)])])):hs("",!0)])]),E("p",Kv,X(T.desc),1),E("div",Xv,[E("span",Yv,X(r(T.url)),1),y[2]||(y[2]=E("i",{class:"fas fa-arrow-up-right-from-square res-card__arrow"},null,-1))])],8,Vv))),128))])]))),128))])])])]))}}),Jv=Cn(Qv,[["__scopeId","data-v-fdc07a4f"]]),Zv="/assets/avatar-DQvqWlfS.png",s1={class:"page-section about-view"},n1={class:"about-head"},a1={class:"about-head__identity"},e1={class:"about-head__avatar"},t1=["src"],l1={key:1,class:"about-head__initial"},o1={class:"about-head__names"},p1={class:"about-head__name"},c1={class:"about-head__role"},r1={key:0,class:"about-intro"},i1={class:"about-body"},u1={key:0,class:"about-section"},h1={class:"about-section__title"},d1={class:"timeline"},m1={class:"tl-year num"},j1={class:"tl-content"},f1={class:"tl-title"},g1={key:0,class:"tl-desc"},b1={key:1,class:"about-grid"},_1={class:"about-cell__title"},y1={key:0,class:"about-cell__item-name"},v1={key:1,class:"about-cell__item-desc"},w1={class:"about-foot"},C1={class:"about-foot__contacts"},k1=["href"],x1={class:"about-foot__thanks"},P1=Gs({name:"AboutView",__name:"About",setup(s){const{t:n,locale:a}=nn(),e=os({}),t=Object.assign({"/src/assets/avatar.png":Zv}),l=ss(()=>Object.values(t)[0]||""),o=ss(()=>n("thanks")),p=un.author,c=p.trim().charAt(0).toUpperCase(),u=ss(()=>e.value?.introduction||""),r=ss(()=>e.value?.experience||[]),i=ss(()=>e.value?.section||[]);function d(){try{e.value=n_()||{}}catch(h){console.error("Failed to load about data:",h),e.value={introduction:"",experience:[],section:[],contacts:[]}}}return d(),Hs(a,()=>{d()}),(h,y)=>(N(),B("div",s1,[E("header",n1,[E("div",a1,[E("div",e1,[l.value?(N(),B("img",{key:0,src:l.value,alt:"avatar",class:"about-head__avatar-img",draggable:"false"},null,8,t1)):(N(),B("span",l1,X(cs(c)),1))]),E("div",o1,[E("h1",p1,X(cs(p)),1),E("p",c1,X(cs(n)("developer")),1)])]),u.value?(N(),B("p",r1,X(u.value),1)):hs("",!0)]),E("main",i1,[r.value.length?(N(),B("section",u1,[E("div",h1,X(cs(n)("experience")),1),E("div",d1,[(N(!0),B(bs,null,Ms(r.value,(j,T)=>(N(),B("div",{key:T,class:"tl-item"},[E("div",m1,X(j.year),1),E("div",j1,[E("div",f1,X(j.title),1),j.desc?(N(),B("div",g1,X(j.desc),1)):hs("",!0)])]))),128))])])):hs("",!0),i.value.length?(N(),B("div",b1,[(N(!0),B(bs,null,Ms(i.value,(j,T)=>(N(),B("div",{key:T,class:"about-cell"},[E("div",_1,X(j.title),1),(N(!0),B(bs,null,Ms(j.items,(g,m)=>(N(),B("div",{key:m,class:"about-cell__item"},[g.name?(N(),B("span",y1,X(g.name),1)):hs("",!0),g.desc?(N(),B("span",v1,X(g.desc),1)):hs("",!0)]))),128))]))),128))])):hs("",!0)]),E("footer",w1,[E("div",C1,[(N(!0),B(bs,null,Ms(e.value.contacts,j=>(N(),B("a",{key:j.label,href:j.link,target:"_blank",rel:"noopener noreferrer",class:"about-foot__contact"},[E("i",{class:Us(j.icon)},null,2),E("span",null,X(j.value),1)],8,k1))),128))]),E("p",x1,X(o.value),1)])]))}}),E1=Cn(P1,[["__scopeId","data-v-4ad5ac55"]]),S1=["innerHTML"],gc=`
  <svg width="16" height="16" viewBox="0 0 14 14" fill="currentColor">
    <path d="M3 2C2.44772 2 2 2.44772 2 3V9C2 9.55228 2.44772 10 3 10H9C9.55228 10 10 9.55228 10 9V3C10 2.44772 9.55228 2 9 2H3ZM1 3C1 1.89543 1.89543 1 3 1H9C10.1046 1 11 1.89543 11 3V9C11 10.1046 10.1046 11 9 11H3C1.89543 11 1 10.1046 1 9V3Z"/>
    <path d="M5 4C4.44772 4 4 4.44772 4 5V11C4 11.5523 4.44772 12 5 12H11C11.5523 12 12 11.5523 12 11V5C12 4.44772 11.5523 4 11 4H5Z"/>
  </svg>
`,T1=`
  <svg width="16" height="16" viewBox="0 0 14 14" fill="currentColor">
    <path d="M11.3536 3.64645C11.5488 3.84171 11.5488 4.15829 11.3536 4.35355L5.35355 10.3536C5.15829 10.5488 4.84171 10.5488 4.64645 10.3536L2.64645 8.35355C2.45118 8.15829 2.45118 7.84171 2.64645 7.64645C2.84171 7.45118 3.15829 7.45118 3.35355 7.64645L5 9.29289L10.6464 3.64645C10.8417 3.45118 11.1583 3.45118 11.3536 3.64645Z"/>
  </svg>
`,A1=Gs({__name:"RenderMarkdown",props:{rawMarkdown:{default:""},articlePath:{default:""},articleTitle:{default:""}},emits:["markdown-rendered"],setup(s,{emit:n}){const{t:a}=nn(),e=Object.assign({}),t=s,l=n,o=os(""),p=Je("markdownContainer");async function c(m){const w=u(m,t.articlePath);o.value=w,await Js(),l("markdown-rendered"),h(),y(),r()}function u(m,w){try{const f=w.replace(/^[./]*/,"").replace(/\.md$/,"").split("/").slice(0,-1).join("/"),x=k=>{if(/^(https?:)?\/\//i.test(k)||k.startsWith("/"))return k;const _=(f+"/"+k).split("/").filter(J=>J&&J!=="."),D=[];_.forEach(J=>J===".."?D.pop():D.push(J));const S=D.join("/"),I=[`@data/content-src/${S}`,`${S}`];for(const J of I){const U=Object.keys(e).find(F=>F.endsWith(`/${S}`)||F===J);if(U)return e[U]}return k};return m.replace(/<img\s+([^>]*?)src=["']([^"']+)["'](.*?)>/gi,(k,_,D,S)=>`<img ${_}src="${x(D.trim())}"${S}>`)}catch(f){return console.warn("rewriteImageLinks failed",f),m}}function r(){const m=p.value;m&&(i(m),d(m))}function i(m){if(!t.articleTitle)return;const w=t.articleTitle.trim().toLowerCase();m.querySelectorAll("h1").forEach(f=>{f.textContent.trim().toLowerCase()===w?f.remove():f.replaceWith(Object.assign(document.createElement("h2"),{...Object.fromEntries(Array.from(f.attributes).map(k=>[k.name,k.value])),innerHTML:f.innerHTML}))})}function d(m){const w=(f,x)=>{const k=window.scrollY+f.getBoundingClientRect().top-88;window.scrollTo({top:Math.max(0,k),behavior:"smooth"}),setTimeout(()=>x.blur(),300)};m.querySelectorAll("h2, h3, h4, h5, h6").forEach(f=>{f.querySelector(".heading-anchor")?.remove();const x=Object.assign(document.createElement("button"),{type:"button",className:"heading-anchor",textContent:"#",ariaLabel:a("anchorHeading"),tabIndex:0,ariaHidden:"false"});x.addEventListener("click",k=>{k.stopPropagation(),w(f,x)}),x.addEventListener("keydown",k=>{(k.key==="Enter"||k.key===" ")&&(k.preventDefault(),w(f,x))}),f.appendChild(x)})}function h(){const m=p.value;m&&m.querySelectorAll("pre").forEach(w=>{if(w.querySelector(".code-block-header"))return;const f=w.querySelector("code");if(!f)return;const x=(f.className.match(/language-(\w+)/)||["","text"])[1],k=document.createElement("div");k.className="code-block-header d-flex align-items-center justify-content-between";const _=document.createElement("span");_.className="code-language",_.textContent=x;const D=document.createElement("button");D.type="button",D.className="copy-button btn-icon d-flex align-items-center justify-content-center",D.setAttribute("aria-label",a("copyCode")),D.innerHTML=gc,D.addEventListener("click",()=>T(f.textContent??"",D)),k.append(_,D);const S=document.createElement("div");S.className="code-block-wrapper",w.parentNode?.insertBefore(S,w),S.append(k,w)})}function y(){const m=p.value;m&&m.querySelectorAll("table").forEach(w=>{if(w.closest(".table-copyable"))return;const f=document.createElement("div");f.className="table-copyable";const x=document.createElement("button");x.type="button",x.className="table-copy-btn",x.setAttribute("aria-label",a("copyTable")),x.innerHTML=gc,x.addEventListener("click",()=>j(w,x)),f.append(x),w.parentNode?.insertBefore(f,w),f.append(w)})}function j(m,w){const x=Array.from(m.querySelectorAll("tr")).map(k=>Array.from(k.querySelectorAll("th, td")).map(_=>(_.textContent||"").trim().replace(/\s+/g," ")).join("	"));T(x.join(`
`),w)}async function T(m,w){try{await navigator.clipboard.writeText(m)}catch(f){console.error(a("copyFailed"),f);const x=document.createElement("textarea");x.value=m,document.body.appendChild(x),x.select(),document.execCommand("copy"),document.body.removeChild(x)}finally{g(w)}}function g(m){const w=m.innerHTML;m.style.color="var(--primary)",m.innerHTML=T1,setTimeout(()=>{m.innerHTML=w,m.style.color=""},1200)}return Hs(()=>t.rawMarkdown,m=>c(m),{immediate:!0}),(m,w)=>(N(),B("div",{class:"markdown-body",innerHTML:o.value,ref_key:"markdownContainer",ref:p},null,8,S1))}}),R1={class:"on-this-page"},L1={class:"otp-header"},M1={class:"otp-title"},D1={class:"otp-list"},O1=["onClick","onKeydown"],I1={class:"otp-text"},N1={key:0,class:"otp-sublist"},F1=["onClick","onKeydown"],$1={class:"otp-subtext"},B1=Gs({__name:"OnThisPage",props:{containerSelector:{default:".markdown-body"},levels:{default:()=>[2,3,4,5,6]},offset:{default:8}},emits:["navigate"],setup(s,{expose:n,emit:a}){const e=s,t=a,{t:l}=nn(),o=os([]),p=os(""),c=os(null),u=os(null),r=os(null),i=os(null);function d(){c.value&&(c.value.disconnect(),c.value=null),u.value&&(clearTimeout(u.value),u.value=null),r.value&&(clearInterval(r.value),r.value=null),i.value&&(i.value.disconnect(),i.value=null)}function h(){o.value=[],p.value="",d(),Js(()=>j())}function y(){g(),m()}function j(){setTimeout(()=>y(),100);const f=()=>{const x=document.querySelector(e.containerSelector);x&&(r.value&&(clearInterval(r.value),r.value=null),c.value||(c.value=new MutationObserver(()=>{clearTimeout(u.value??void 0),u.value=window.setTimeout(()=>y(),100)}),c.value.observe(x,{childList:!0,subtree:!0,attributes:!0,attributeFilter:["id"]})),y())};f(),!r.value&&!document.querySelector(e.containerSelector)&&(r.value=window.setInterval(f,200))}function T(f){try{const x=f.cloneNode(!0);return x.querySelectorAll(".heading-anchor")?.forEach(k=>k.remove()),(x.textContent||"").replace(/\s*#\s*$/,"").trim()}catch{return(f.textContent||"").replace(/\s*#\s*$/,"").trim()}}function g(){const f=document.querySelector(e.containerSelector);if(!f){o.value=[];return}const x=e.levels.map(U=>`h${U}`).join(","),k=Array.from(f.querySelectorAll(x));if(k.length===0){o.value=[];return}k.forEach(U=>{if(!U.id){let F=U.textContent.trim().toLowerCase().replace(/[^\u4e00-\u9fa5a-zA-Z0-9\s-]/g,"").replace(/\s+/g,"-");F||(F=`section-${Math.random().toString(36).substring(2,9)}`);let G=F,is=1;for(;document.getElementById(G);)G=`${F}-${is++}`;U.id=G}});const _=new Set(e.levels),D=Math.min(...e.levels),S=D+1,I=[];let J=null;for(const U of k){const F=parseInt(U.tagName.substring(1),10);if(!_.has(F))continue;const G={id:U.id,text:T(U),level:F,children:[]};F===D?(I.push(G),J=G):J&&F>=S?J.children.push(G):I.push(G)}o.value=I}function m(){i.value&&(i.value.disconnect(),i.value=null);const f=document.querySelector(e.containerSelector);if(!f)return;const x=e.levels.map(_=>`h${_}`).join(","),k=Array.from(f.querySelectorAll(x));k.length!==0&&(i.value=new IntersectionObserver(_=>{for(const D of _)D.isIntersecting&&(p.value=D.target.id)},{rootMargin:`-${e.offset}px 0px -60% 0px`,threshold:0}),k.forEach(_=>i.value?.observe(_)))}function w(f){t("navigate",f);const x=document.getElementById(f);if(!x)return;const k=x.getBoundingClientRect().top+window.scrollY-e.offset,_=typeof window.matchMedia=="function"&&window.matchMedia("(prefers-reduced-motion: reduce)").matches,D=()=>{window.scrollTo({top:k,behavior:_?"auto":"smooth"})};try{(document.body&&document.body.style&&document.body.style.overflow)==="hidden"?setTimeout(D,80):D()}catch{D()}}return wn(()=>{g(),m(),Js(()=>{j()})}),qn(()=>{d()}),n({refreshToc:y,resetToc:h}),(f,x)=>(N(),B("nav",R1,[E("div",L1,[E("span",M1,X(cs(l)("tableOfContents")),1)]),E("ul",D1,[(N(!0),B(bs,null,Ms(o.value,(k,_)=>(N(),B("li",{key:_,class:Us(["otp-item",{active:p.value===k.id}])},[E("a",{class:"otp-link",role:"button",tabindex:"0",onClick:_n(D=>w(k.id),["prevent"]),onKeydown:fp(_n(D=>w(k.id),["prevent"]),["enter"])},[E("span",I1,X(k.text),1)],40,O1),k.children&&k.children.length?(N(),B("ul",N1,[(N(!0),B(bs,null,Ms(k.children,(D,S)=>(N(),B("li",{key:S,class:Us(["otp-subitem",{active:p.value===D.id}])},[E("a",{class:"otp-sublink",role:"button",tabindex:"0",onClick:_n(I=>w(D.id),["prevent"]),onKeydown:fp(_n(I=>w(D.id),["prevent"]),["enter"])},[E("span",$1,X(D.text),1)],40,F1)],2))),128))])):hs("",!0)],2))),128))])]))}}),Qi=Cn(B1,[["__scopeId","data-v-309aab29"]]),q1=["aria-label"],z1={class:"offcanvas-panel offcanvas-right border-start rounded-0 shadow-sm"},V1={class:"offcanvas-section"},U1={class:"offcanvas-card"},bc="toc",H1=Gs({__name:"TocDrawer",setup(s){const{t:n}=nn(),a=os(!1),e=os(null),t=os(!1),l=os(!1),o=os(0),p=os(0),c=os((typeof window<"u"?window.innerHeight:1024)-160),u=os(!1),r=os(!1),i=os(null);function d(){return{gap:48,minTop:20,maxTop:Math.max(0,window.innerHeight-40-20-48)}}function h(S){const{minTop:I,maxTop:J}=d();return Math.max(I,Math.min(J,S))}function y(S){e.value=S,!a.value&&(a.value=!0,requestAnimationFrame(()=>{window.dispatchEvent(new CustomEvent("floating-buttons-base-top",{detail:{baseTop:e.value,source:bc}})),a.value=!1}))}function j(S){const I=S.detail,J=I?.baseTop;if(I?.source===bc)return;const{gap:F}=d();if(typeof J=="number"){const G=J-F;c.value=h(G)}}function T(){const S=window.innerWidth<992;t.value=S,c.value=h(c.value)}function g(){l.value||(_(),r.value=!0)}function m(){r.value=!1,D()}function w(S){S.preventDefault(),l.value=!0,u.value=!1,o.value=S.touches[0].clientY,p.value=c.value}function f(S){u.value=!0;const J=S.touches[0].clientY-o.value,U=h(p.value+J);c.value=U;const{gap:F}=d();y(U+F),S.preventDefault()}function x(S){S.preventDefault(),l.value=!1,u.value||g()}function k(){m()}function _(){try{const S=window.scrollY||window.pageYOffset||document.documentElement.scrollTop||0;i.value=S;const I=document.body;I&&(I.style.position="fixed",I.style.top=`-${S}px`,I.style.left="0",I.style.right="0",I.style.overflow="hidden"),document.documentElement&&(document.documentElement.style.overscrollBehavior="contain")}catch{const S=document.documentElement,I=document.body;S&&(S.style.overflow="hidden",S.style.overscrollBehavior="contain"),I&&(I.style.overflow="hidden",I.style.overscrollBehavior="contain")}}function D(){try{const S=document.body;S&&(S.style.position="",S.style.top="",S.style.left="",S.style.right="",S.style.overflow=""),document.documentElement&&(document.documentElement.style.overscrollBehavior=""),typeof i.value=="number"&&(window.scrollTo(0,i.value),i.value=null)}catch{const S=document.documentElement,I=document.body;S&&(S.style.overflow="",S.style.overscrollBehavior=""),I&&(I.style.overflow="",I.style.overscrollBehavior="")}}return wn(()=>{window.addEventListener("resize",T,{passive:!0}),window.addEventListener("floating-buttons-base-top",j),T(),y(c.value+48)}),qn(()=>{window.removeEventListener("resize",T),window.removeEventListener("floating-buttons-base-top",j),D()}),(S,I)=>(N(),B(bs,null,[Ua(E("button",{class:"toc-drawer-btn d-lg-none d-flex align-items-center justify-content-center",onClick:g,"aria-label":cs(n)("openToc"),onTouchstart:_n(w,["prevent","stop"]),onTouchmove:_n(f,["prevent","stop"]),onTouchend:_n(x,["prevent","stop"]),style:ba({top:c.value+"px"})},[...I[0]||(I[0]=[E("i",{class:"fas fa-bookmark"},null,-1)])],44,q1),[[qr,t.value]]),r.value?(N(),B("div",{key:0,class:"mobile-offcanvas d-lg-none",onClick:_n(m,["self"])},[E("div",z1,[E("div",V1,[E("div",U1,[ws(Qi,{containerSelector:".markdown-body",levels:[2,3],offset:88,onNavigate:k})])])]),E("div",{class:"offcanvas-backdrop",onClick:m})])):hs("",!0)],64))}}),G1=Cn(H1,[["__scopeId","data-v-3f1455a9"]]),W1={class:"container view-container article-view"},K1={class:"row py-4 px-0"},X1={class:"col-12 col-lg-2 order-2 order-lg-1 docs-sidebar-col"},Y1={key:0,class:"navigation-container mb-0"},Q1={class:"col-12 col-lg-8 order-1 order-lg-2 docs-main-col",ref:"mainContent"},J1={class:"article-panel"},Z1={class:"article-panel__body"},sw={class:"article-content"},nw={key:0,class:"article-meta"},aw={class:"article-title"},ew={class:"article-meta__row"},tw={key:0,class:"meta-line"},lw={key:1,class:"meta-line"},ow=["aria-label"],pw={key:0,class:"article-meta__tags"},cw=["onClick"],rw={key:2,class:"article-navigation"},iw={key:0,class:"article-nav-spacer","aria-hidden":"true"},uw={class:"nav-details"},hw={class:"nav-label"},dw={class:"nav-title"},mw={class:"nav-details"},jw={class:"nav-label"},fw={class:"nav-title"},gw={class:"col-12 col-lg-2 order-3 docs-toc-col"},bw={key:0,class:"toc-container mt-0"},_w=Gs({name:"ArticleView",head(){return{title:this.currentPost?`${this.currentPost.title} - zorrooz’s blog`:"zorrooz’s blog",meta:this.currentPost?.description?[{name:"description",content:this.currentPost.description}]:[]}},__name:"Article",props:{path:{type:[String,Array],default:""}},setup(s){const{t:n,locale:a}=nn(),e=Ka(),t=Ie(),l=os(""),o=os(""),p=os([]),c=os([]),u=os(typeof window<"u"?window.innerWidth:1024),r=os(0),i=os(!1);let d=null,h=!1;const y=Je("onThisPageRef"),j=Je("leftSidebarContent"),T=Je("rightSidebarContent"),g=ss(()=>u.value>=992),m=ss(()=>!!(k.value?.path&&k.value.path.startsWith("notes/"))),w=ss(()=>n("updatedAt")),f=ss(()=>n("prevPage")),x=ss(()=>n("nextPage")),k=ss(()=>o.value?p.value.find(ts=>{const rs=ts.path.replace(/\.md$/,""),fs=o.value.replace(/\.md$/,"");return!!(rs===fs||fs.endsWith("-en")&&rs===fs.replace(/-en$/,"")||!fs.endsWith("-en")&&rs===fs+"-en")}):null),_=ss(()=>{if(!k.value)return[];const[ts,rs]=k.value.path.replace(/\.md$/,"").split("/"),fs=[],_s=(O,W)=>{if(!W)return;const H=String(W).replace(/^\/+/,"").split("/"),Q=H[0]==="article"?1:0,gs=H[Q],C=H[Q+1];if(gs!==ts||C!==rs)return;const P=H.slice(Q+2),R=`${gs}/${C}/${P.join("/")}`;fs.push({title:O,path:`${R}.md`})};if(Array.isArray(c.value)){for(const O of c.value)if(Array.isArray(O.items))for(const W of O.items)W?.name===rs&&(Array.isArray(W.articles)&&W.articles.forEach(H=>_s(H.title,H.articleUrl)),Array.isArray(W.categories)&&W.categories.forEach(H=>{Array.isArray(H.articles)&&H.articles.forEach(Q=>_s(Q.title,Q.articleUrl))}))}return fs}),D=ss(()=>k.value?_.value.findIndex(ts=>ts.path.replace(/\.md$/,"")===k.value.path.replace(/\.md$/,"")):-1),S=ss(()=>{const ts=D.value;return ts>0?_.value[ts-1]:null}),I=ss(()=>{const ts=D.value,rs=_.value.length-1;return ts>=0&&ts<rs?_.value[ts+1]:null}),J=ss(()=>{const ts=k.value?.wordCount;return typeof ts=="number"?Yi(ts):0});function U(ts){return n("articleReadingTime",{minutes:ts})}function F(ts){return Qn(`/article/${ts.replace(/\.md$/,"")}`)}function G(ts){if(!ts)return;const rs=a.value==="en-US"?"/en":"/zh";t.push({path:`${rs}/`,query:{tag:ts,from:e.fullPath}})}async function is(){if(!k.value)return;let ts="";try{const rs=await e_(k.value.path);rs&&(ts=rs.trim())}catch{}if(!ts){const rs=document.querySelector(".markdown-body");if(!rs)return;const fs=rs.cloneNode(!0);fs.querySelectorAll(".heading-anchor, .code-block-header, .copy-button, .table-copy-btn").forEach(O=>O.remove());const _s=(fs.innerText||"").trim();ts=`${k.value.title}

${_s}`}try{await navigator.clipboard.writeText(ts)}catch(rs){console.error(n("copyFailed"),rs);const fs=document.createElement("textarea");fs.value=ts,document.body.appendChild(fs),fs.select(),document.execCommand("copy"),document.body.removeChild(fs)}finally{i.value=!0,d&&clearTimeout(d),d=window.setTimeout(()=>{i.value=!1},1200)}}function ps(){const ts=[],rs=(fs,_s,O=[],W="",H=0)=>{if(typeof _s!="string"||!_s.trim())return;const Q=_s.replace(/^\/+/,"").split("/"),gs=Q[0]==="article"?1:0,C=Q[gs],P=Q[gs+1],R=Q.slice(gs+2),$=R[0]||"__root__",V=R.join("/"),q=`${C}/${P}/${V}`,b={title:fs,path:`${q}.md`,date:W||"",tags:Array.isArray(O)?O:[],preview:"",category:`${C}/${P}/${$}`,wordCount:typeof H=="number"?H:0};ts.push(b)};try{const fs=Ft();Array.isArray(fs)&&fs.forEach(_s=>{Array.isArray(_s.items)&&_s.items.forEach(O=>{const W=O?.stats?.latestDate||"";Array.isArray(O.articles)&&O.articles.forEach(H=>rs(H.title,H.articleUrl,H?.tags||[],W,H?.wordCount)),Array.isArray(O.categories)&&O.categories.forEach(H=>{const Q=H?.stats?.latestDate||W||"";Array.isArray(H.articles)&&H.articles.forEach(gs=>rs(gs.title,gs.articleUrl,gs?.tags||[],Q,gs?.wordCount))})})}),p.value=ts,c.value=fs}catch{p.value=[],c.value=[]}}function Y(){u.value=window.innerWidth,Rs()}function us(){h||(h=!0,requestAnimationFrame(()=>{h=!1;const rs=document.documentElement.scrollHeight-window.innerHeight;r.value=rs>0?Math.min(100,window.scrollY/rs*100):0}))}function Ns(ts){return Array.isArray(ts)?ts.join("/"):typeof ts=="string"&&ts.length?ts:""}function Xs(){l.value="";try{const ts=Ns(e.params.path),rs=a.value==="en-US";let fs=p.value.find(_s=>{const O=_s.path.replace(/\.md$/,"");return rs?O===ts||O===ts.replace(/-en$/,"")+"-en":O===ts||O===ts.replace(/-en$/,"")});if(fs||(fs=p.value.find(_s=>{const O=_s.path.replace(/\.md$/,""),W=ts.replace(/-en$/,"");return O===W||O===W+"-en"})),!fs)throw new Error(`Article not found: ${ts}`);o.value=fs.path,l.value=Qb(o.value),Js(()=>{typeof window>"u"||(window.scrollTo({top:0,behavior:"smooth"}),Rs(),y.value?.refreshToc())})}catch{l.value=`# Article Not Found

The requested article could not be loaded. Please check the URL.`,Js(()=>{typeof window>"u"||(window.scrollTo({top:0,behavior:"smooth"}),y.value?.refreshToc())})}}function Rs(){if(typeof window>"u")return;const rs=document.querySelector("header")?.offsetHeight||60,fs=window.innerHeight,_s=Math.max(200,fs-rs-24-24);[j.value,T.value].forEach(W=>{W&&(W.style.maxHeight=`${_s}px`,W.style.overflowY="auto")})}function zs(){Rs(),Js(()=>y.value?.refreshToc())}return ps(),Xs(),Hs(a,(ts,rs)=>{ts!==rs&&(ps(),Xs())}),Hs(()=>e.params.path,(ts,rs)=>{const fs=Ns(rs),_s=Ns(ts);fs!==_s&&(y.value?.resetToc(),Xs())}),Hs(l,()=>{Js(()=>Rs())}),wn(()=>{u.value=window.innerWidth,window.addEventListener("resize",Y),window.addEventListener("scroll",us,{passive:!0})}),qn(()=>{window.removeEventListener("resize",Y),window.removeEventListener("scroll",us),d&&clearTimeout(d)}),(ts,rs)=>{const fs=Le("router-link");return N(),B("div",W1,[E("div",{class:"reading-progress",style:ba({width:r.value+"%"}),"aria-hidden":"true"},null,4),E("div",K1,[E("div",X1,[E("div",{class:"sticky-sidebar",ref_key:"leftSidebarContent",ref:j},[g.value?(N(),B("div",Y1,[ws(Xi)])):hs("",!0)],512)]),E("div",Q1,[E("div",J1,[E("div",Z1,[E("div",sw,[k.value?(N(),B("div",nw,[E("h1",aw,X(k.value.title),1),E("div",ew,[m.value&&k.value.date?(N(),B("span",tw,[rs[0]||(rs[0]=E("i",{class:"fas fa-calendar-alt"},null,-1)),As(X(w.value)+" "+X(k.value.date),1)])):hs("",!0),J.value>0?(N(),B("span",lw,[rs[1]||(rs[1]=E("i",{class:"fas fa-clock"},null,-1)),As(X(U(J.value)),1)])):hs("",!0),E("button",{type:"button",class:Us(["article-copy-btn",{"article-copy-btn--copied":i.value}]),onClick:is,"aria-label":cs(n)("copyArticle"),"aria-live":"polite"},[E("i",{class:Us(i.value?"fas fa-check":"fas fa-copy")},null,2),E("span",null,X(i.value?cs(n)("copied"):cs(n)("copyArticle")),1)],10,ow)]),k.value.tags?.length?(N(),B("div",pw,[(N(!0),B(bs,null,Ms(k.value.tags,(_s,O)=>(N(),B("span",{key:O,class:"article-tag",onClick:W=>G(_s)}," # "+X(_s),9,cw))),128))])):hs("",!0)])):hs("",!0),l.value?(N(),Sa(A1,{key:1,rawMarkdown:l.value,articlePath:k.value?.path||"",articleTitle:k.value?.title||"",onMarkdownRendered:zs},null,8,["rawMarkdown","articlePath","articleTitle"])):hs("",!0),l.value?(N(),B("nav",rw,[!S.value&&I.value?(N(),B("span",iw)):hs("",!0),S.value?(N(),Sa(fs,{key:1,to:F(S.value.path),class:"article-nav-item prev"},{default:vn(()=>[rs[2]||(rs[2]=E("div",{class:"nav-arrow"},[E("i",{class:"fas fa-arrow-left"})],-1)),E("div",uw,[E("div",hw,X(f.value),1),E("div",dw,X(S.value.title),1)])]),_:1},8,["to"])):hs("",!0),I.value?(N(),Sa(fs,{key:2,to:F(I.value.path),class:"article-nav-item next"},{default:vn(()=>[E("div",mw,[E("div",jw,X(x.value),1),E("div",fw,X(I.value.title),1)]),rs[3]||(rs[3]=E("div",{class:"nav-arrow"},[E("i",{class:"fas fa-arrow-right"})],-1))]),_:1},8,["to"])):hs("",!0)])):hs("",!0)])])])],512),E("div",gw,[E("div",{class:"sticky-sidebar",ref_key:"rightSidebarContent",ref:T},[g.value?(N(),B("div",bw,[ws(Qi,{ref_key:"onThisPageRef",ref:y,containerSelector:".markdown-body",levels:[2,3],offset:88},null,512)])):hs("",!0)],512)])]),l.value?(N(),Sa(G1,{key:0})):hs("",!0)])}}}),yw=Cn(_w,[["__scopeId","data-v-bc45a872"]]),ae=()=>{if(typeof window<"u"){const s=localStorage.getItem("locale");if(s===Sn[1])return ua[1];if(s===Sn[0])return ua[0];if(navigator.language&&navigator.language.toLowerCase().startsWith("en"))return ua[1]}return ua[0]},_c=s=>[{path:`${s}/`,name:`${s}-Home`,component:av},{path:`${s}/category`,name:`${s}-Category`,component:Sv},{path:`${s}/resource`,name:`${s}-Resource`,component:Jv},{path:`${s}/about`,name:`${s}-About`,component:E1},{path:`${s}/article/:path*`,name:`${s}-Article`,component:yw,props:!0}],vw=[{path:"/",redirect:()=>ae()},{path:"/category",redirect:s=>`${ae()}${s.path}`},{path:"/resource",redirect:s=>`${ae()}${s.path}`},{path:"/about",redirect:s=>`${ae()}${s.path}`},{path:"/article/:path*",redirect:s=>`${ae()}${s.path}`},..._c(ua[0]),..._c(ua[1])],ww={mounted(s){if(typeof window>"u"||window.matchMedia("(prefers-reduced-motion: reduce)").matches||!("IntersectionObserver"in window))return;s.classList.add("reveal");const n=new IntersectionObserver(([a])=>{a.isIntersecting&&(s.classList.add("reveal-visible"),n.disconnect())},{threshold:.12});n.observe(s),s.__revealIo=n},unmounted(s){s.__revealIo?.disconnect(),s.__revealIo=null}};ys(()=>import("./bootstrap.esm-inoyDHN5.js"),[]);rj(ry,{routes:vw},({app:s,router:n,initialState:a,isClient:e,routePath:t})=>{const l=ij();s.use(l),a.pinia&&(l.state.value=a.pinia),s.use(n),s.use(_t),s.mixin(nm),s.directive("reveal",ww);const o=wo();o.initTheme(),e&&o.initLocale(),e&&"serviceWorker"in navigator&&window.addEventListener("load",()=>{navigator.serviceWorker.register("/sw.js").catch(()=>{})})});export{bs as F,Cw as T,ys as _,Ie as a,E as b,Sa as c,Gs as d,cs as e,Ua as f,B as g,As as h,hs as i,Ms as j,fp as k,_n as l,Qn as m,Js as n,N as o,Cn as p,os as r,X as t,nn as u,kw as v,Hs as w};
