const __vite__mapDeps=(i,m=__vite__mapDeps,d=(m.f||(m.f=["assets/SearchModal-BJTHk7sV.js","assets/SearchModal-CcvjA6bn.css"])))=>i.map(i=>d[i]);
(function(){const n=document.createElement("link").relList;if(n&&n.supports&&n.supports("modulepreload"))return;for(const t of document.querySelectorAll('link[rel="modulepreload"]'))e(t);new MutationObserver(t=>{for(const l of t)if(l.type==="childList")for(const o of l.addedNodes)o.tagName==="LINK"&&o.rel==="modulepreload"&&e(o)}).observe(document,{childList:!0,subtree:!0});function a(t){const l={};return t.integrity&&(l.integrity=t.integrity),t.referrerPolicy&&(l.referrerPolicy=t.referrerPolicy),t.crossOrigin==="use-credentials"?l.credentials="include":t.crossOrigin==="anonymous"?l.credentials="omit":l.credentials="same-origin",l}function e(t){if(t.ep)return;t.ep=!0;const l=a(t);fetch(t.href,l)}})();const Gr="modulepreload",Ur=function(s){return"/"+s},ql={},Mp=function(n,a,e){let t=Promise.resolve();if(a&&a.length>0){let c=function(u){return Promise.all(u.map(r=>Promise.resolve(r).then(i=>({status:"fulfilled",value:i}),i=>({status:"rejected",reason:i}))))};document.getElementsByTagName("link");const o=document.querySelector("meta[property=csp-nonce]"),p=o?.nonce||o?.getAttribute("nonce");t=c(a.map(u=>{if(u=Ur(u),u in ql)return;ql[u]=!0;const r=u.endsWith(".css"),i=r?'[rel="stylesheet"]':"";if(document.querySelector(`link[href="${u}"]${i}`))return;const d=document.createElement("link");if(d.rel=r?"stylesheet":Gr,r||(d.as="script"),d.crossOrigin="",d.href=u,p&&d.setAttribute("nonce",p),document.head.appendChild(d),r)return new Promise((h,j)=>{d.addEventListener("load",h),d.addEventListener("error",()=>j(new Error(`Unable to preload CSS for ${u}`)))})}))}function l(o){const p=new Event("vite:preloadError",{cancelable:!0});if(p.payload=o,window.dispatchEvent(p),!p.defaultPrevented)throw o}return t.then(o=>{for(const p of o||[])p.status==="rejected"&&l(p.reason);return n().catch(l)})};function Mt(s,n={},a){for(const e in s){const t=s[e],l=a?`${a}:${e}`:e;typeof t=="object"&&t!==null?Mt(t,n,l):typeof t=="function"&&(n[l]=t)}return n}const Dp=(()=>{if(console.createTask)return console.createTask;const s={run:n=>n()};return()=>s})();function Lp(s,n,a,e){for(let t=a;t<s.length;t+=1)try{const l=e?e.run(()=>s[t](...n)):s[t](...n);if(l&&typeof l.then=="function")return Promise.resolve(l).then(()=>Lp(s,n,t+1,e))}catch(l){return Promise.reject(l)}}function Wr(s,n,a){if(s.length>0)return Lp(s,n,0,Dp(a))}function Kr(s,n,a){if(s.length>0){const e=Dp(a);return Promise.all(s.map(t=>e.run(()=>t(...n))))}}function ht(s,n){for(const a of[...s])a(n)}var Yr=class{_hooks;_before;_after;_deprecatedHooks;_deprecatedMessages;constructor(){this._hooks={},this._before=void 0,this._after=void 0,this._deprecatedMessages=void 0,this._deprecatedHooks={},this.hook=this.hook.bind(this),this.callHook=this.callHook.bind(this),this.callHookWith=this.callHookWith.bind(this)}hook(s,n,a={}){if(!s||typeof n!="function")return()=>{};const e=s;let t;for(;this._deprecatedHooks[s];)t=this._deprecatedHooks[s],s=t.to;if(t&&!a.allowDeprecated){let l=t.message;l||(l=`${e} hook has been deprecated`+(t.to?`, please use ${t.to}`:"")),this._deprecatedMessages||(this._deprecatedMessages=new Set),this._deprecatedMessages.has(l)||(console.warn(l),this._deprecatedMessages.add(l))}if(!n.name)try{Object.defineProperty(n,"name",{get:()=>"_"+s.replace(/\W+/g,"_")+"_hook_cb",configurable:!0})}catch{}return this._hooks[s]=this._hooks[s]||[],this._hooks[s].push(n),()=>{n&&(this.removeHook(s,n),n=void 0)}}hookOnce(s,n){let a,e=(...t)=>(typeof a=="function"&&a(),a=void 0,e=void 0,n(...t));return a=this.hook(s,e),a}removeHook(s,n){const a=this._hooks[s];if(a){const e=a.indexOf(n);e!==-1&&a.splice(e,1),a.length===0&&(this._hooks[s]=void 0)}}clearHook(s){this._hooks[s]=void 0}deprecateHook(s,n){this._deprecatedHooks[s]=typeof n=="string"?{to:n}:n;const a=this._hooks[s]||[];this._hooks[s]=void 0;for(const e of a)this.hook(s,e)}deprecateHooks(s){for(const n in s)this.deprecateHook(n,s[n])}addHooks(s){const n=Mt(s),a=Object.keys(n).map(e=>this.hook(e,n[e]));return()=>{for(const e of a)e();a.length=0}}removeHooks(s){const n=Mt(s);for(const a in n)this.removeHook(a,n[a])}removeAllHooks(){this._hooks={}}callHook(s,...n){return this.callHookWith(Wr,s,n)}callHookParallel(s,...n){return this.callHookWith(Kr,s,n)}callHookWith(s,n,a){const e=this._before||this._after?{name:n,args:a,context:{}}:void 0;this._before&&ht(this._before,e);const t=s(this._hooks[n]?[...this._hooks[n]]:[],a,n);return t instanceof Promise?t.finally(()=>{this._after&&e&&ht(this._after,e)}):(this._after&&e&&ht(this._after,e),t)}beforeEach(s){return this._before=this._before||[],this._before.push(s),()=>{if(this._before!==void 0){const n=this._before.indexOf(s);n!==-1&&this._before.splice(n,1)}}}afterEach(s){return this._after=this._after||[],this._after.push(s),()=>{if(this._after!==void 0){const n=this._after.indexOf(s);n!==-1&&this._after.splice(n,1)}}}};function Xr(){return new Yr}const Jr=new Set(["link","style","script","noscript"]),Qr=new Set(["title","titleTemplate","script","style","noscript"]),Dt=new Set(["base","meta","link","style","script","noscript"]),Zr=new Set(["title","base","htmlAttrs","bodyAttrs","meta","link","style","script","noscript"]),si=new Set(["base","title","titleTemplate","bodyAttrs","htmlAttrs","templateParams"]),ni=new Set(["key","tagPosition","tagPriority","tagDuplicateStrategy","innerHTML","textContent","processTemplateParams"]),ai=new Set(["templateParams","htmlAttrs","bodyAttrs"]),ei=new Set(["theme-color","google-site-verification","og","article","book","profile","twitter","author"]),ti=["name","property","http-equiv"],li=new Set(["viewport","description","keywords","robots"]);function Op(s){const n=s.split(":");return n.length?ei.has(n[1]):!1}function Lt(s){const{props:n,tag:a}=s;if(si.has(a))return a;if(a==="link"&&n.rel==="canonical")return"canonical";if(a==="link"&&n.rel==="alternate"){if(n.hreflang)return`alternate:${n.hreflang}`;if(n.type)return`alternate:${n.type}:${n.href||""}`}if(n.charset)return"charset";if(s.tag==="meta"){for(const e of ti)if(n[e]!==void 0){const t=n[e],l=t&&typeof t=="string"&&t.includes(":"),o=t&&li.has(t),c=!(l||o)&&s.key?`:key:${s.key}`:"";return`${a}:${t}${c}`}}if(s.key)return`${a}:key:${s.key}`;if(n.id)return`${a}:id:${n.id}`;if(a==="link"&&n.rel==="alternate")return`alternate:${n.href||""}`;if(Qr.has(a)){const e=s.textContent||s.innerHTML;if(e)return`${a}:content:${e}`}}function Fp(s){const n=s._h||s._d;if(n)return n;const a=s.textContent||s.innerHTML;return a||`${s.tag}:${Object.entries(s.props).map(([e,t])=>`${e}:${String(t)}`).join(",")}`}function Oe(s,n,a){typeof s==="function"&&(!a||a!=="titleTemplate"&&!(a[0]==="o"&&a[1]==="n"))&&(s=s());const t=n?n(a,s):s;if(Array.isArray(t))return t.map(l=>Oe(l,n));if(t?.constructor===Object){const l={};for(const o of Object.keys(t))l[o]=Oe(t[o],n,o);return l}return t}function oi(s,n){const a=s==="style"?new Map:new Set;function e(t){if(t==null||t===void 0)return;const l=String(t).trim();if(l)if(s==="style"){const[o,...p]=l.split(":").map(c=>c?c.trim():"");o&&p.length&&a.set(o,p.join(":"))}else l.split(" ").filter(Boolean).forEach(o=>a.add(o))}return typeof n=="string"?s==="style"?n.split(";").forEach(e):e(n):Array.isArray(n)?n.forEach(t=>e(t)):n&&typeof n=="object"&&Object.entries(n).forEach(([t,l])=>{l&&l!=="false"&&(s==="style"?a.set(String(t).trim(),String(l)):e(t))}),a}function $p(s,n){if(s.props=s.props||{},!n)return s;if(s.tag==="templateParams")return s.props=n,s;const a=Dt.has(s.tag)||s.tag==="htmlAttrs"||s.tag==="bodyAttrs";return Object.entries(n).forEach(([e,t])=>{if(e==="__proto__"||e==="constructor"||e==="prototype")return;if(t===null){s.props[e]=null;return}if(e==="class"||e==="style"){s.props[e]=oi(e,t);return}if(ni.has(e)){if((e==="textContent"||e==="innerHTML")&&typeof t=="object"){let u=n.type;if(n.type||(u="application/json"),!u?.endsWith("json")&&u!=="speculationrules")return;n.type=u,s.props.type=u,s[e]=JSON.stringify(t)}else s[e]=t;return}const l=e.startsWith("data-"),o=a&&!l?e.toLowerCase():e,p=String(t),c=s.tag==="meta"&&o==="content";p==="true"||p===""?s.props[o]=l||c?p:!0:!t&&l&&p==="false"?s.props[o]="false":t!==void 0&&(s.props[o]=t)}),s}function pi(s,n){const a=typeof n=="object"&&typeof n!="function"?n:{[s==="script"||s==="noscript"||s==="style"?"innerHTML":"textContent"]:n},e=$p({tag:s,props:{}},a);return e.key&&Jr.has(e.tag)&&(e.props["data-hid"]=e._h=e.key),e.tag==="script"&&typeof e.innerHTML=="object"&&(e.innerHTML=JSON.stringify(e.innerHTML),e.props.type=e.props.type||"application/json"),Array.isArray(e.props.content)?e.props.content.map(t=>({...e,props:{...e.props,content:t}})):e}function ci(s,n){if(!s)return[];typeof s=="function"&&(s=s());const a=(t,l)=>{for(let o=0;o<n.length;o++)l=n[o](t,l);return l};s=a(void 0,s);const e=[];return s=Oe(s,a),Object.entries(s||{}).forEach(([t,l])=>{if(l!==void 0)for(const o of Array.isArray(l)?l:[l])e.push(pi(t,o))}),e.flat()}const Hl=(s,n)=>s._w===n._w?s._p-n._p:s._w-n._w,Vl={base:-10,title:10},ri={critical:-8,high:-1,low:2},Gl={meta:{"content-security-policy":-30,charset:-20,viewport:-15},link:{preconnect:20,stylesheet:60,preload:70,modulepreload:70,prefetch:90,"dns-prefetch":90,prerender:90},script:{async:30,defer:80,sync:50},style:{imported:40,sync:60}},ii=/@import/,La=s=>s===""||s===!0;function ui(s,n){if(typeof n.tagPriority=="number")return n.tagPriority;let a=100;const e=ri[n.tagPriority]||0,t=s.resolvedOptions.disableCapoSorting?{link:{},script:{},style:{}}:Gl;if(n.tag in Vl)a=Vl[n.tag];else if(n.tag==="meta"){const l=n.props["http-equiv"]==="content-security-policy"?"content-security-policy":n.props.charset?"charset":n.props.name==="viewport"?"viewport":null;l&&(a=Gl.meta[l])}else if(n.tag==="link"&&n.props.rel)a=t.link[n.props.rel];else if(n.tag==="script"){const l=String(n.props.type);La(n.props.async)?a=t.script.async:n.props.src&&!La(n.props.defer)&&!La(n.props.async)&&l!=="module"&&!l.endsWith("json")||n.innerHTML&&!l.endsWith("json")?a=t.script.sync:(La(n.props.defer)&&n.props.src&&!La(n.props.async)||l==="module")&&(a=t.script.defer)}else n.tag==="style"&&(a=n.innerHTML&&ii.test(n.innerHTML)?t.style.imported:t.style.sync);return(a||100)+e}function Ul(s,n){const a=typeof n=="function"?n(s):n,e=a.key||String(s.plugins.size+1);s.plugins.get(e)||(s.plugins.set(e,a),s.hooks.addHooks(a.hooks||{}))}function hi(s={}){const n=Xr();n.addHooks(s.hooks||{});const a=!s.document,e=new Map,t=new Map,l=new Set,o={_entryCount:1,plugins:t,dirty:!1,resolvedOptions:s,hooks:n,ssr:a,entries:e,headEntries(){return[...e.values()]},use:p=>Ul(o,p),push(p,c){const u={...c||{}};delete u.head;const r=u._index??o._entryCount++,i={_i:r,input:p,options:u},d={_poll(h=!1){o.dirty=!0,!h&&l.add(r),n.callHook("entries:updated",o)},dispose(){e.delete(r)&&o.invalidate()},patch(h){(!u.mode||u.mode==="server"&&a||u.mode==="client"&&!a)&&(i.input=h,e.set(r,i),d._poll())}};return d.patch(p),d},async resolveTags(){const p={tagMap:new Map,tags:[],entries:[...o.entries.values()]};for(await n.callHook("entries:resolve",p);l.size;){const d=l.values().next().value;l.delete(d);const h=e.get(d);if(h){const j={tags:ci(h.input,s.propResolvers||[]).map(m=>Object.assign(m,h.options)),entry:h};await n.callHook("entries:normalize",j),h._tags=j.tags.map((m,P)=>(m._w=ui(o,m),m._p=(h._i<<10)+P,m._d=Lt(m),m._d||(m._h=Fp(m)),m))}}let c=!1;p.entries.flatMap(d=>(d._tags||[]).map(h=>({...h,props:{...h.props}}))).sort(Hl).reduce((d,h)=>{const j=h._d||h._h;if(!d.has(j))return d.set(j,h);const m=d.get(j);if((h?.tagDuplicateStrategy||(ai.has(h.tag)?"merge":null)||(h.key&&h.key===m.key?"merge":null))==="merge"){const C={...m.props};Object.entries(h.props).forEach(([b,x])=>C[b]=b==="style"?new Map([...m.props.style||new Map,...x]):b==="class"?new Set([...m.props.class||new Set,...x]):x),d.set(j,{...h,props:C})}else h._p>>10===m._p>>10&&h.tag==="meta"&&Op(j)?(d.set(j,Object.assign([...Array.isArray(m)?m:[m],h],h)),c=!0):(h._w===m._w?h._p>m._p:h?._w<m?._w)&&d.set(j,h);return d},p.tagMap);const u=p.tagMap.get("title"),r=p.tagMap.get("titleTemplate");if(o._title=u?.textContent,r){const d=r?.textContent;if(o._titleTemplate=d,d){let h=typeof d=="function"?d(u?.textContent):d;typeof h=="string"&&!o.plugins.has("template-params")&&(h=h.replace("%s",u?.textContent||"")),u?h===null?p.tagMap.delete("title"):p.tagMap.set("title",{...u,textContent:h}):(r.tag="title",r.textContent=h)}}p.tags=Array.from(p.tagMap.values()),c&&(p.tags=p.tags.flat().sort(Hl)),await n.callHook("tags:beforeResolve",p),await n.callHook("tags:resolve",p),await n.callHook("tags:afterResolve",p);const i=[];for(const d of p.tags){const{innerHTML:h,tag:j,props:m}=d;if(Zr.has(j)&&!(Object.keys(m).length===0&&!d.innerHTML&&!d.textContent)&&!(j==="meta"&&!m.content&&!m["http-equiv"]&&!m.charset)){if(j==="script"&&h){if(String(m.type).endsWith("json")){const P=typeof h=="string"?h:JSON.stringify(h);d.innerHTML=P.replace(/</g,"\\u003C")}else typeof h=="string"&&(d.innerHTML=h.replace(new RegExp(`</${j}`,"g"),`<\\/${j}`));d._d=Lt(d)}i.push(d)}}return i},invalidate(){for(const p of e.values())l.add(p._i);o.dirty=!0,n.callHook("entries:updated",o)}};return(s?.plugins||[]).forEach(p=>Ul(o,p)),o.hooks.callHook("init",o),s.init?.forEach(p=>p&&o.push(p)),o}async function Ip(s,n={}){const a=n.document||s.resolvedOptions.document;if(!a||!s.dirty)return;const e={shouldRender:!0,tags:[]};if(await s.hooks.callHook("dom:beforeRender",e),!!e.shouldRender)return s._domUpdatePromise||(s._domUpdatePromise=new Promise(async t=>{const l=new Map,o=new Promise(h=>{s.resolveTags().then(j=>{h(j.map(m=>{const P=l.get(m._d)||0,C={tag:m,id:(P?`${m._d}:${P}`:m._d)||m._h,shouldRender:!0};return m._d&&Op(m._d)&&l.set(m._d,P+1),C}))})});let p=s._dom;if(!p){p={title:a.title,elMap:new Map().set("htmlAttrs",a.documentElement).set("bodyAttrs",a.body)};for(const h of["body","head"]){const j=a[h]?.children;for(const m of j){const P=m.tagName.toLowerCase();if(!Dt.has(P))continue;const C=$p({tag:P,props:{}},{innerHTML:m.innerHTML,...m.getAttributeNames().reduce((b,x)=>(b[x]=m.getAttribute(x),b),{})||{}});if(C.key=m.getAttribute("data-hid")||void 0,C._d=Lt(C)||Fp(C),p.elMap.has(C._d)){let b=1,x=C._d;for(;p.elMap.has(x);)x=`${C._d}:${b++}`;p.elMap.set(x,m)}else p.elMap.set(C._d,m)}}}p.pendingSideEffects={...p.sideEffects},p.sideEffects={};function c(h,j,m){const P=`${h}:${j}`;p.sideEffects[P]=m,delete p.pendingSideEffects[P]}function u({id:h,$el:j,tag:m}){const P=m.tag.endsWith("Attrs");p.elMap.set(h,j),P||(m.textContent&&m.textContent!==j.textContent&&(j.textContent=m.textContent),m.innerHTML&&m.innerHTML!==j.innerHTML&&(j.innerHTML=m.innerHTML),c(h,"el",()=>{j?.remove(),p.elMap.delete(h)}));for(const C in m.props){if(!Object.prototype.hasOwnProperty.call(m.props,C))continue;const b=m.props[C];if(C.startsWith("on")&&typeof b=="function"){const f=j?.dataset;if(f&&f[`${C}fired`]){const k=C.slice(0,-5);b.call(j,new Event(k.substring(2)))}j.getAttribute(`data-${C}`)!==""&&((m.tag==="bodyAttrs"?a.defaultView:j).addEventListener(C.substring(2),b.bind(j)),j.setAttribute(`data-${C}`,""));continue}const x=`attr:${C}`;if(C==="class"){if(!b)continue;for(const f of b)P&&c(h,`${x}:${f}`,()=>j.classList.remove(f)),!j.classList.contains(f)&&j.classList.add(f)}else if(C==="style"){if(!b)continue;for(const[f,k]of b)c(h,`${x}:${f}`,()=>{j.style.removeProperty(f)}),j.style.setProperty(f,k)}else b!==!1&&b!==null&&(j.getAttribute(C)!==b&&j.setAttribute(C,b===!0?"":String(b)),P&&c(h,x,()=>j.removeAttribute(C)))}}const r=[],i={bodyClose:void 0,bodyOpen:void 0,head:void 0},d=await o;for(const h of d){const{tag:j,shouldRender:m,id:P}=h;if(m){if(j.tag==="title"){a.title=j.textContent,c("title","",()=>a.title=p.title);continue}h.$el=h.$el||p.elMap.get(P),h.$el?u(h):Dt.has(j.tag)&&r.push(h)}}for(const h of r){const j=h.tag.tagPosition||"head";h.$el=a.createElement(h.tag.tag),u(h),i[j]=i[j]||a.createDocumentFragment(),i[j].appendChild(h.$el)}for(const h of d)await s.hooks.callHook("dom:renderTag",h,a,c);i.head&&a.head.appendChild(i.head),i.bodyOpen&&a.body.insertBefore(i.bodyOpen,a.body.firstChild),i.bodyClose&&a.body.appendChild(i.bodyClose);for(const h in p.pendingSideEffects)p.pendingSideEffects[h]();s._dom=p,await s.hooks.callHook("dom:rendered",{renders:d}),t()}).finally(()=>{s._domUpdatePromise=void 0,s.dirty=!1})),s._domUpdatePromise}function di(s={}){const n=s.domOptions?.render||Ip;s.document=s.document||(typeof window<"u"?document:void 0);const a=s.document?.head.querySelector('script[id="unhead:payload"]')?.innerHTML||!1;return hi({...s,plugins:[...s.plugins||[],{key:"client",hooks:{"entries:updated":n}}],init:[a?JSON.parse(a):!1,...s.init||[]]})}function mi(s,n){let a=0;return()=>{const e=++a;n(()=>{a===e&&s()})}}/**
* @vue/shared v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/function ll(s){const n=Object.create(null);for(const a of s.split(","))n[a]=1;return a=>a in n}const ks={},Ca=[],An=()=>{},Np=()=>!1,Xe=s=>s.charCodeAt(0)===111&&s.charCodeAt(1)===110&&(s.charCodeAt(2)>122||s.charCodeAt(2)<97),ol=s=>s.startsWith("onUpdate:"),Gs=Object.assign,pl=(s,n)=>{const a=s.indexOf(n);a>-1&&s.splice(a,1)},ji=Object.prototype.hasOwnProperty,xs=(s,n)=>ji.call(s,n),is=Array.isArray,ka=s=>Je(s)==="[object Map]",Bp=s=>Je(s)==="[object Set]",ds=s=>typeof s=="function",Ms=s=>typeof s=="string",ta=s=>typeof s=="symbol",Rs=s=>s!==null&&typeof s=="object",zp=s=>(Rs(s)||ds(s))&&ds(s.then)&&ds(s.catch),qp=Object.prototype.toString,Je=s=>qp.call(s),gi=s=>Je(s).slice(8,-1),Hp=s=>Je(s)==="[object Object]",cl=s=>Ms(s)&&s!=="NaN"&&s[0]!=="-"&&""+parseInt(s,10)===s,qa=ll(",key,ref,ref_for,ref_key,onVnodeBeforeMount,onVnodeMounted,onVnodeBeforeUpdate,onVnodeUpdated,onVnodeBeforeUnmount,onVnodeUnmounted"),Qe=s=>{const n=Object.create(null);return(a=>n[a]||(n[a]=s(a)))},fi=/-\w/g,yn=Qe(s=>s.replace(fi,n=>n.slice(1).toUpperCase())),bi=/\B([A-Z])/g,la=Qe(s=>s.replace(bi,"-$1").toLowerCase()),Ze=Qe(s=>s.charAt(0).toUpperCase()+s.slice(1)),dt=Qe(s=>s?`on${Ze(s)}`:""),Zn=(s,n)=>!Object.is(s,n),Pe=(s,...n)=>{for(let a=0;a<s.length;a++)s[a](...n)},Vp=(s,n,a,e=!1)=>{Object.defineProperty(s,n,{configurable:!0,enumerable:!1,writable:e,value:a})},rl=s=>{const n=parseFloat(s);return isNaN(n)?s:n},_i=s=>{const n=Ms(s)?Number(s):NaN;return isNaN(n)?s:n};let Wl;const st=()=>Wl||(Wl=typeof globalThis<"u"?globalThis:typeof self<"u"?self:typeof window<"u"?window:typeof global<"u"?global:{});function oa(s){if(is(s)){const n={};for(let a=0;a<s.length;a++){const e=s[a],t=Ms(e)?xi(e):oa(e);if(t)for(const l in t)n[l]=t[l]}return n}else if(Ms(s)||Rs(s))return s}const yi=/;(?![^(]*\))/g,vi=/:([^]+)/,wi=/\/\*[^]*?\*\//g;function xi(s){const n={};return s.replace(wi,"").split(yi).forEach(a=>{if(a){const e=a.split(vi);e.length>1&&(n[e[0].trim()]=e[1].trim())}}),n}function Ns(s){let n="";if(Ms(s))n=s;else if(is(s))for(let a=0;a<s.length;a++){const e=Ns(s[a]);e&&(n+=e+" ")}else if(Rs(s))for(const a in s)s[a]&&(n+=a+" ");return n.trim()}const Ci="itemscope,allowfullscreen,formnovalidate,ismap,nomodule,novalidate,readonly",ki=ll(Ci);function Gp(s){return!!s||s===""}const Up=s=>!!(s&&s.__v_isRef===!0),G=s=>Ms(s)?s:s==null?"":is(s)||Rs(s)&&(s.toString===qp||!ds(s.toString))?Up(s)?G(s.value):JSON.stringify(s,Wp,2):String(s),Wp=(s,n)=>Up(n)?Wp(s,n.value):ka(n)?{[`Map(${n.size})`]:[...n.entries()].reduce((a,[e,t],l)=>(a[mt(e,l)+" =>"]=t,a),{})}:Bp(n)?{[`Set(${n.size})`]:[...n.values()].map(a=>mt(a))}:ta(n)?mt(n):Rs(n)&&!is(n)&&!Hp(n)?String(n):n,mt=(s,n="")=>{var a;return ta(s)?`Symbol(${(a=s.description)!=null?a:n})`:s};/**
* @vue/reactivity v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/let Js;class Kp{constructor(n=!1){this.detached=n,this._active=!0,this._on=0,this.effects=[],this.cleanups=[],this._isPaused=!1,this.parent=Js,!n&&Js&&(this.index=(Js.scopes||(Js.scopes=[])).push(this)-1)}get active(){return this._active}pause(){if(this._active){this._isPaused=!0;let n,a;if(this.scopes)for(n=0,a=this.scopes.length;n<a;n++)this.scopes[n].pause();for(n=0,a=this.effects.length;n<a;n++)this.effects[n].pause()}}resume(){if(this._active&&this._isPaused){this._isPaused=!1;let n,a;if(this.scopes)for(n=0,a=this.scopes.length;n<a;n++)this.scopes[n].resume();for(n=0,a=this.effects.length;n<a;n++)this.effects[n].resume()}}run(n){if(this._active){const a=Js;try{return Js=this,n()}finally{Js=a}}}on(){++this._on===1&&(this.prevScope=Js,Js=this)}off(){this._on>0&&--this._on===0&&(Js=this.prevScope,this.prevScope=void 0)}stop(n){if(this._active){this._active=!1;let a,e;for(a=0,e=this.effects.length;a<e;a++)this.effects[a].stop();for(this.effects.length=0,a=0,e=this.cleanups.length;a<e;a++)this.cleanups[a]();if(this.cleanups.length=0,this.scopes){for(a=0,e=this.scopes.length;a<e;a++)this.scopes[a].stop(!0);this.scopes.length=0}if(!this.detached&&this.parent&&!n){const t=this.parent.scopes.pop();t&&t!==this&&(this.parent.scopes[this.index]=t,t.index=this.index)}this.parent=void 0}}}function il(s){return new Kp(s)}function Yp(){return Js}function Si(s,n=!1){Js&&Js.cleanups.push(s)}let As;const jt=new WeakSet;class Xp{constructor(n){this.fn=n,this.deps=void 0,this.depsTail=void 0,this.flags=5,this.next=void 0,this.cleanup=void 0,this.scheduler=void 0,Js&&Js.active&&Js.effects.push(this)}pause(){this.flags|=64}resume(){this.flags&64&&(this.flags&=-65,jt.has(this)&&(jt.delete(this),this.trigger()))}notify(){this.flags&2&&!(this.flags&32)||this.flags&8||Qp(this)}run(){if(!(this.flags&1))return this.fn();this.flags|=2,Kl(this),Zp(this);const n=As,a=vn;As=this,vn=!0;try{return this.fn()}finally{sc(this),As=n,vn=a,this.flags&=-3}}stop(){if(this.flags&1){for(let n=this.deps;n;n=n.nextDep)dl(n);this.deps=this.depsTail=void 0,Kl(this),this.onStop&&this.onStop(),this.flags&=-2}}trigger(){this.flags&64?jt.add(this):this.scheduler?this.scheduler():this.runIfDirty()}runIfDirty(){Ot(this)&&this.run()}get dirty(){return Ot(this)}}let Jp=0,Ha,Va;function Qp(s,n=!1){if(s.flags|=8,n){s.next=Va,Va=s;return}s.next=Ha,Ha=s}function ul(){Jp++}function hl(){if(--Jp>0)return;if(Va){let n=Va;for(Va=void 0;n;){const a=n.next;n.next=void 0,n.flags&=-9,n=a}}let s;for(;Ha;){let n=Ha;for(Ha=void 0;n;){const a=n.next;if(n.next=void 0,n.flags&=-9,n.flags&1)try{n.trigger()}catch(e){s||(s=e)}n=a}}if(s)throw s}function Zp(s){for(let n=s.deps;n;n=n.nextDep)n.version=-1,n.prevActiveLink=n.dep.activeLink,n.dep.activeLink=n}function sc(s){let n,a=s.depsTail,e=a;for(;e;){const t=e.prevDep;e.version===-1?(e===a&&(a=t),dl(e),Pi(e)):n=e,e.dep.activeLink=e.prevActiveLink,e.prevActiveLink=void 0,e=t}s.deps=n,s.depsTail=a}function Ot(s){for(let n=s.deps;n;n=n.nextDep)if(n.dep.version!==n.version||n.dep.computed&&(nc(n.dep.computed)||n.dep.version!==n.version))return!0;return!!s._dirty}function nc(s){if(s.flags&4&&!(s.flags&16)||(s.flags&=-17,s.globalVersion===Za)||(s.globalVersion=Za,!s.isSSR&&s.flags&128&&(!s.deps&&!s._dirty||!Ot(s))))return;s.flags|=2;const n=s.dep,a=As,e=vn;As=s,vn=!0;try{Zp(s);const t=s.fn(s._value);(n.version===0||Zn(t,s._value))&&(s.flags|=128,s._value=t,n.version++)}catch(t){throw n.version++,t}finally{As=a,vn=e,sc(s),s.flags&=-3}}function dl(s,n=!1){const{dep:a,prevSub:e,nextSub:t}=s;if(e&&(e.nextSub=t,s.prevSub=void 0),t&&(t.prevSub=e,s.nextSub=void 0),a.subs===s&&(a.subs=e,!e&&a.computed)){a.computed.flags&=-5;for(let l=a.computed.deps;l;l=l.nextDep)dl(l,!0)}!n&&!--a.sc&&a.map&&a.map.delete(a.key)}function Pi(s){const{prevDep:n,nextDep:a}=s;n&&(n.nextDep=a,s.prevDep=void 0),a&&(a.prevDep=n,s.nextDep=void 0)}let vn=!0;const ac=[];function zn(){ac.push(vn),vn=!1}function qn(){const s=ac.pop();vn=s===void 0?!0:s}function Kl(s){const{cleanup:n}=s;if(s.cleanup=void 0,n){const a=As;As=void 0;try{n()}finally{As=a}}}let Za=0;class Ai{constructor(n,a){this.sub=n,this.dep=a,this.version=a.version,this.nextDep=this.prevDep=this.nextSub=this.prevSub=this.prevActiveLink=void 0}}class ml{constructor(n){this.computed=n,this.version=0,this.activeLink=void 0,this.subs=void 0,this.map=void 0,this.key=void 0,this.sc=0,this.__v_skip=!0}track(n){if(!As||!vn||As===this.computed)return;let a=this.activeLink;if(a===void 0||a.sub!==As)a=this.activeLink=new Ai(As,this),As.deps?(a.prevDep=As.depsTail,As.depsTail.nextDep=a,As.depsTail=a):As.deps=As.depsTail=a,ec(a);else if(a.version===-1&&(a.version=this.version,a.nextDep)){const e=a.nextDep;e.prevDep=a.prevDep,a.prevDep&&(a.prevDep.nextDep=e),a.prevDep=As.depsTail,a.nextDep=void 0,As.depsTail.nextDep=a,As.depsTail=a,As.deps===a&&(As.deps=e)}return a}trigger(n){this.version++,Za++,this.notify(n)}notify(n){ul();try{for(let a=this.subs;a;a=a.prevSub)a.sub.notify()&&a.sub.dep.notify()}finally{hl()}}}function ec(s){if(s.dep.sc++,s.sub.flags&4){const n=s.dep.computed;if(n&&!s.dep.subs){n.flags|=20;for(let e=n.deps;e;e=e.nextDep)ec(e)}const a=s.dep.subs;a!==s&&(s.prevSub=a,a&&(a.nextSub=s)),s.dep.subs=s}}const Fe=new WeakMap,ja=Symbol(""),Ft=Symbol(""),se=Symbol("");function Qs(s,n,a){if(vn&&As){let e=Fe.get(s);e||Fe.set(s,e=new Map);let t=e.get(a);t||(e.set(a,t=new ml),t.map=e,t.key=a),t.track()}}function On(s,n,a,e,t,l){const o=Fe.get(s);if(!o){Za++;return}const p=c=>{c&&c.trigger()};if(ul(),n==="clear")o.forEach(p);else{const c=is(s),u=c&&cl(a);if(c&&a==="length"){const r=Number(e);o.forEach((i,d)=>{(d==="length"||d===se||!ta(d)&&d>=r)&&p(i)})}else switch((a!==void 0||o.has(void 0))&&p(o.get(a)),u&&p(o.get(se)),n){case"add":c?u&&p(o.get("length")):(p(o.get(ja)),ka(s)&&p(o.get(Ft)));break;case"delete":c||(p(o.get(ja)),ka(s)&&p(o.get(Ft)));break;case"set":ka(s)&&p(o.get(ja));break}}hl()}function Ti(s,n){const a=Fe.get(s);return a&&a.get(n)}function _a(s){const n=ys(s);return n===s?n:(Qs(n,"iterate",se),bn(s)?n:n.map(Ws))}function nt(s){return Qs(s=ys(s),"iterate",se),s}const Ri={__proto__:null,[Symbol.iterator](){return gt(this,Symbol.iterator,Ws)},concat(...s){return _a(this).concat(...s.map(n=>is(n)?_a(n):n))},entries(){return gt(this,"entries",s=>(s[1]=Ws(s[1]),s))},every(s,n){return Rn(this,"every",s,n,void 0,arguments)},filter(s,n){return Rn(this,"filter",s,n,a=>a.map(Ws),arguments)},find(s,n){return Rn(this,"find",s,n,Ws,arguments)},findIndex(s,n){return Rn(this,"findIndex",s,n,void 0,arguments)},findLast(s,n){return Rn(this,"findLast",s,n,Ws,arguments)},findLastIndex(s,n){return Rn(this,"findLastIndex",s,n,void 0,arguments)},forEach(s,n){return Rn(this,"forEach",s,n,void 0,arguments)},includes(...s){return ft(this,"includes",s)},indexOf(...s){return ft(this,"indexOf",s)},join(s){return _a(this).join(s)},lastIndexOf(...s){return ft(this,"lastIndexOf",s)},map(s,n){return Rn(this,"map",s,n,void 0,arguments)},pop(){return Oa(this,"pop")},push(...s){return Oa(this,"push",s)},reduce(s,...n){return Yl(this,"reduce",s,n)},reduceRight(s,...n){return Yl(this,"reduceRight",s,n)},shift(){return Oa(this,"shift")},some(s,n){return Rn(this,"some",s,n,void 0,arguments)},splice(...s){return Oa(this,"splice",s)},toReversed(){return _a(this).toReversed()},toSorted(s){return _a(this).toSorted(s)},toSpliced(...s){return _a(this).toSpliced(...s)},unshift(...s){return Oa(this,"unshift",s)},values(){return gt(this,"values",Ws)}};function gt(s,n,a){const e=nt(s),t=e[n]();return e!==s&&!bn(s)&&(t._next=t.next,t.next=()=>{const l=t._next();return l.done||(l.value=a(l.value)),l}),t}const Ei=Array.prototype;function Rn(s,n,a,e,t,l){const o=nt(s),p=o!==s&&!bn(s),c=o[n];if(c!==Ei[n]){const i=c.apply(s,l);return p?Ws(i):i}let u=a;o!==s&&(p?u=function(i,d){return a.call(this,Ws(i),d,s)}:a.length>2&&(u=function(i,d){return a.call(this,i,d,s)}));const r=c.call(o,u,e);return p&&t?t(r):r}function Yl(s,n,a,e){const t=nt(s);let l=a;return t!==s&&(bn(s)?a.length>3&&(l=function(o,p,c){return a.call(this,o,p,c,s)}):l=function(o,p,c){return a.call(this,o,Ws(p),c,s)}),t[n](l,...e)}function ft(s,n,a){const e=ys(s);Qs(e,"iterate",se);const t=e[n](...a);return(t===-1||t===!1)&&fl(a[0])?(a[0]=ys(a[0]),e[n](...a)):t}function Oa(s,n,a=[]){zn(),ul();const e=ys(s)[n].apply(s,a);return hl(),qn(),e}const Mi=ll("__proto__,__v_isRef,__isVue"),tc=new Set(Object.getOwnPropertyNames(Symbol).filter(s=>s!=="arguments"&&s!=="caller").map(s=>Symbol[s]).filter(ta));function Di(s){ta(s)||(s=String(s));const n=ys(this);return Qs(n,"has",s),n.hasOwnProperty(s)}class lc{constructor(n=!1,a=!1){this._isReadonly=n,this._isShallow=a}get(n,a,e){if(a==="__v_skip")return n.__v_skip;const t=this._isReadonly,l=this._isShallow;if(a==="__v_isReactive")return!t;if(a==="__v_isReadonly")return t;if(a==="__v_isShallow")return l;if(a==="__v_raw")return e===(t?l?Hi:rc:l?cc:pc).get(n)||Object.getPrototypeOf(n)===Object.getPrototypeOf(e)?n:void 0;const o=is(n);if(!t){let c;if(o&&(c=Ri[a]))return c;if(a==="hasOwnProperty")return Di}const p=Reflect.get(n,a,Es(n)?n:e);if((ta(a)?tc.has(a):Mi(a))||(t||Qs(n,"get",a),l))return p;if(Es(p)){const c=o&&cl(a)?p:p.value;return t&&Rs(c)?It(c):c}return Rs(p)?t?It(p):re(p):p}}class oc extends lc{constructor(n=!1){super(!1,n)}set(n,a,e,t){let l=n[a];if(!this._isShallow){const c=na(l);if(!bn(e)&&!na(e)&&(l=ys(l),e=ys(e)),!is(n)&&Es(l)&&!Es(e))return c||(l.value=e),!0}const o=is(n)&&cl(a)?Number(a)<n.length:xs(n,a),p=Reflect.set(n,a,e,Es(n)?n:t);return n===ys(t)&&(o?Zn(e,l)&&On(n,"set",a,e):On(n,"add",a,e)),p}deleteProperty(n,a){const e=xs(n,a);n[a];const t=Reflect.deleteProperty(n,a);return t&&e&&On(n,"delete",a,void 0),t}has(n,a){const e=Reflect.has(n,a);return(!ta(a)||!tc.has(a))&&Qs(n,"has",a),e}ownKeys(n){return Qs(n,"iterate",is(n)?"length":ja),Reflect.ownKeys(n)}}class Li extends lc{constructor(n=!1){super(!0,n)}set(n,a){return!0}deleteProperty(n,a){return!0}}const Oi=new oc,Fi=new Li,$i=new oc(!0);const $t=s=>s,be=s=>Reflect.getPrototypeOf(s);function Ii(s,n,a){return function(...e){const t=this.__v_raw,l=ys(t),o=ka(l),p=s==="entries"||s===Symbol.iterator&&o,c=s==="keys"&&o,u=t[s](...e),r=a?$t:n?$e:Ws;return!n&&Qs(l,"iterate",c?Ft:ja),{next(){const{value:i,done:d}=u.next();return d?{value:i,done:d}:{value:p?[r(i[0]),r(i[1])]:r(i),done:d}},[Symbol.iterator](){return this}}}}function _e(s){return function(...n){return s==="delete"?!1:s==="clear"?void 0:this}}function Ni(s,n){const a={get(t){const l=this.__v_raw,o=ys(l),p=ys(t);s||(Zn(t,p)&&Qs(o,"get",t),Qs(o,"get",p));const{has:c}=be(o),u=n?$t:s?$e:Ws;if(c.call(o,t))return u(l.get(t));if(c.call(o,p))return u(l.get(p));l!==o&&l.get(t)},get size(){const t=this.__v_raw;return!s&&Qs(ys(t),"iterate",ja),t.size},has(t){const l=this.__v_raw,o=ys(l),p=ys(t);return s||(Zn(t,p)&&Qs(o,"has",t),Qs(o,"has",p)),t===p?l.has(t):l.has(t)||l.has(p)},forEach(t,l){const o=this,p=o.__v_raw,c=ys(p),u=n?$t:s?$e:Ws;return!s&&Qs(c,"iterate",ja),p.forEach((r,i)=>t.call(l,u(r),u(i),o))}};return Gs(a,s?{add:_e("add"),set:_e("set"),delete:_e("delete"),clear:_e("clear")}:{add(t){!n&&!bn(t)&&!na(t)&&(t=ys(t));const l=ys(this);return be(l).has.call(l,t)||(l.add(t),On(l,"add",t,t)),this},set(t,l){!n&&!bn(l)&&!na(l)&&(l=ys(l));const o=ys(this),{has:p,get:c}=be(o);let u=p.call(o,t);u||(t=ys(t),u=p.call(o,t));const r=c.call(o,t);return o.set(t,l),u?Zn(l,r)&&On(o,"set",t,l):On(o,"add",t,l),this},delete(t){const l=ys(this),{has:o,get:p}=be(l);let c=o.call(l,t);c||(t=ys(t),c=o.call(l,t)),p&&p.call(l,t);const u=l.delete(t);return c&&On(l,"delete",t,void 0),u},clear(){const t=ys(this),l=t.size!==0,o=t.clear();return l&&On(t,"clear",void 0,void 0),o}}),["keys","values","entries",Symbol.iterator].forEach(t=>{a[t]=Ii(t,s,n)}),a}function jl(s,n){const a=Ni(s,n);return(e,t,l)=>t==="__v_isReactive"?!s:t==="__v_isReadonly"?s:t==="__v_raw"?e:Reflect.get(xs(a,t)&&t in e?a:e,t,l)}const Bi={get:jl(!1,!1)},zi={get:jl(!1,!0)},qi={get:jl(!0,!1)};const pc=new WeakMap,cc=new WeakMap,rc=new WeakMap,Hi=new WeakMap;function Vi(s){switch(s){case"Object":case"Array":return 1;case"Map":case"Set":case"WeakMap":case"WeakSet":return 2;default:return 0}}function Gi(s){return s.__v_skip||!Object.isExtensible(s)?0:Vi(gi(s))}function re(s){return na(s)?s:gl(s,!1,Oi,Bi,pc)}function ic(s){return gl(s,!1,$i,zi,cc)}function It(s){return gl(s,!0,Fi,qi,rc)}function gl(s,n,a,e,t){if(!Rs(s)||s.__v_raw&&!(n&&s.__v_isReactive))return s;const l=Gi(s);if(l===0)return s;const o=t.get(s);if(o)return o;const p=new Proxy(s,l===2?e:a);return t.set(s,p),p}function sa(s){return na(s)?sa(s.__v_raw):!!(s&&s.__v_isReactive)}function na(s){return!!(s&&s.__v_isReadonly)}function bn(s){return!!(s&&s.__v_isShallow)}function fl(s){return s?!!s.__v_raw:!1}function ys(s){const n=s&&s.__v_raw;return n?ys(n):s}function bl(s){return!xs(s,"__v_skip")&&Object.isExtensible(s)&&Vp(s,"__v_skip",!0),s}const Ws=s=>Rs(s)?re(s):s,$e=s=>Rs(s)?It(s):s;function Es(s){return s?s.__v_isRef===!0:!1}function ns(s){return uc(s,!1)}function _l(s){return uc(s,!0)}function uc(s,n){return Es(s)?s:new Ui(s,n)}class Ui{constructor(n,a){this.dep=new ml,this.__v_isRef=!0,this.__v_isShallow=!1,this._rawValue=a?n:ys(n),this._value=a?n:Ws(n),this.__v_isShallow=a}get value(){return this.dep.track(),this._value}set value(n){const a=this._rawValue,e=this.__v_isShallow||bn(n)||na(n);n=e?n:ys(n),Zn(n,a)&&(this._rawValue=n,this._value=e?n:Ws(n),this.dep.trigger())}}function ls(s){return Es(s)?s.value:s}function Wi(s){return ds(s)?s():ls(s)}const Ki={get:(s,n,a)=>n==="__v_raw"?s:ls(Reflect.get(s,n,a)),set:(s,n,a,e)=>{const t=s[n];return Es(t)&&!Es(a)?(t.value=a,!0):Reflect.set(s,n,a,e)}};function hc(s){return sa(s)?s:new Proxy(s,Ki)}function Yi(s){const n=is(s)?new Array(s.length):{};for(const a in s)n[a]=Ji(s,a);return n}class Xi{constructor(n,a,e){this._object=n,this._key=a,this._defaultValue=e,this.__v_isRef=!0,this._value=void 0}get value(){const n=this._object[this._key];return this._value=n===void 0?this._defaultValue:n}set value(n){this._object[this._key]=n}get dep(){return Ti(ys(this._object),this._key)}}function Ji(s,n,a){const e=s[n];return Es(e)?e:new Xi(s,n,a)}class Qi{constructor(n,a,e){this.fn=n,this.setter=a,this._value=void 0,this.dep=new ml(this),this.__v_isRef=!0,this.deps=void 0,this.depsTail=void 0,this.flags=16,this.globalVersion=Za-1,this.next=void 0,this.effect=this,this.__v_isReadonly=!a,this.isSSR=e}notify(){if(this.flags|=16,!(this.flags&8)&&As!==this)return Qp(this,!0),!0}get value(){const n=this.dep.track();return nc(this),n&&(n.version=this.dep.version),this._value}set value(n){this.setter&&this.setter(n)}}function Zi(s,n,a=!1){let e,t;return ds(s)?e=s:(e=s.get,t=s.set),new Qi(e,t,a)}const ye={},Ie=new WeakMap;let ha;function su(s,n=!1,a=ha){if(a){let e=Ie.get(a);e||Ie.set(a,e=[]),e.push(s)}}function nu(s,n,a=ks){const{immediate:e,deep:t,once:l,scheduler:o,augmentJob:p,call:c}=a,u=f=>t?f:bn(f)||t===!1||t===0?Fn(f,1):Fn(f);let r,i,d,h,j=!1,m=!1;if(Es(s)?(i=()=>s.value,j=bn(s)):sa(s)?(i=()=>u(s),j=!0):is(s)?(m=!0,j=s.some(f=>sa(f)||bn(f)),i=()=>s.map(f=>{if(Es(f))return f.value;if(sa(f))return u(f);if(ds(f))return c?c(f,2):f()})):ds(s)?n?i=c?()=>c(s,2):s:i=()=>{if(d){zn();try{d()}finally{qn()}}const f=ha;ha=r;try{return c?c(s,3,[h]):s(h)}finally{ha=f}}:i=An,n&&t){const f=i,k=t===!0?1/0:t;i=()=>Fn(f(),k)}const P=Yp(),C=()=>{r.stop(),P&&P.active&&pl(P.effects,r)};if(l&&n){const f=n;n=(...k)=>{f(...k),C()}}let b=m?new Array(s.length).fill(ye):ye;const x=f=>{if(!(!(r.flags&1)||!r.dirty&&!f))if(n){const k=r.run();if(t||j||(m?k.some((w,y)=>Zn(w,b[y])):Zn(k,b))){d&&d();const w=ha;ha=r;try{const y=[k,b===ye?void 0:m&&b[0]===ye?[]:b,h];b=k,c?c(n,3,y):n(...y)}finally{ha=w}}}else r.run()};return p&&p(x),r=new Xp(i),r.scheduler=o?()=>o(x,!1):x,h=f=>su(f,!1,r),d=r.onStop=()=>{const f=Ie.get(r);if(f){if(c)c(f,4);else for(const k of f)k();Ie.delete(r)}},n?e?x(!0):b=r.run():o?o(x.bind(null,!0),!0):r.run(),C.pause=r.pause.bind(r),C.resume=r.resume.bind(r),C.stop=C,C}function Fn(s,n=1/0,a){if(n<=0||!Rs(s)||s.__v_skip||(a=a||new Map,(a.get(s)||0)>=n))return s;if(a.set(s,n),n--,Es(s))Fn(s.value,n,a);else if(is(s))for(let e=0;e<s.length;e++)Fn(s[e],n,a);else if(Bp(s)||ka(s))s.forEach(e=>{Fn(e,n,a)});else if(Hp(s)){for(const e in s)Fn(s[e],n,a);for(const e of Object.getOwnPropertySymbols(s))Object.prototype.propertyIsEnumerable.call(s,e)&&Fn(s[e],n,a)}return s}/**
* @vue/runtime-core v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/function ie(s,n,a,e){try{return e?s(...e):s()}catch(t){ue(t,n,a)}}function wn(s,n,a,e){if(ds(s)){const t=ie(s,n,a,e);return t&&zp(t)&&t.catch(l=>{ue(l,n,a)}),t}if(is(s)){const t=[];for(let l=0;l<s.length;l++)t.push(wn(s[l],n,a,e));return t}}function ue(s,n,a,e=!0){const t=n?n.vnode:null,{errorHandler:l,throwUnhandledErrorInProduction:o}=n&&n.appContext.config||ks;if(n){let p=n.parent;const c=n.proxy,u=`https://vuejs.org/error-reference/#runtime-${a}`;for(;p;){const r=p.ec;if(r){for(let i=0;i<r.length;i++)if(r[i](s,c,u)===!1)return}p=p.parent}if(l){zn(),ie(l,null,10,[s,c,u]),qn();return}}au(s,a,t,e,o)}function au(s,n,a,e=!0,t=!1){if(t)throw s;console.error(s)}const tn=[];let Sn=-1;const Sa=[];let Yn=null,va=0;const dc=Promise.resolve();let Ne=null;function Ks(s){const n=Ne||dc;return s?n.then(this?s.bind(this):s):n}function eu(s){let n=Sn+1,a=tn.length;for(;n<a;){const e=n+a>>>1,t=tn[e],l=ne(t);l<s||l===s&&t.flags&2?n=e+1:a=e}return n}function yl(s){if(!(s.flags&1)){const n=ne(s),a=tn[tn.length-1];!a||!(s.flags&2)&&n>=ne(a)?tn.push(s):tn.splice(eu(n),0,s),s.flags|=1,mc()}}function mc(){Ne||(Ne=dc.then(gc))}function tu(s){is(s)?Sa.push(...s):Yn&&s.id===-1?Yn.splice(va+1,0,s):s.flags&1||(Sa.push(s),s.flags|=1),mc()}function Xl(s,n,a=Sn+1){for(;a<tn.length;a++){const e=tn[a];if(e&&e.flags&2){if(s&&e.id!==s.uid)continue;tn.splice(a,1),a--,e.flags&4&&(e.flags&=-2),e(),e.flags&4||(e.flags&=-2)}}}function jc(s){if(Sa.length){const n=[...new Set(Sa)].sort((a,e)=>ne(a)-ne(e));if(Sa.length=0,Yn){Yn.push(...n);return}for(Yn=n,va=0;va<Yn.length;va++){const a=Yn[va];a.flags&4&&(a.flags&=-2),a.flags&8||a(),a.flags&=-2}Yn=null,va=0}}const ne=s=>s.id==null?s.flags&2?-1:1/0:s.id;function gc(s){try{for(Sn=0;Sn<tn.length;Sn++){const n=tn[Sn];n&&!(n.flags&8)&&(n.flags&4&&(n.flags&=-2),ie(n,n.i,n.i?15:14),n.flags&4||(n.flags&=-2))}}finally{for(;Sn<tn.length;Sn++){const n=tn[Sn];n&&(n.flags&=-2)}Sn=-1,tn.length=0,jc(),Ne=null,(tn.length||Sa.length)&&gc()}}let cn=null,fc=null;function Be(s){const n=cn;return cn=s,fc=s&&s.type.__scopeId||null,n}function hn(s,n=cn,a){if(!n||s._n)return s;const e=(...t)=>{e._d&&He(-1);const l=Be(n);let o;try{o=s(...t)}finally{Be(l),e._d&&He(1)}return o};return e._n=!0,e._c=!0,e._d=!0,e}function Aa(s,n){if(cn===null)return s;const a=tt(cn),e=s.dirs||(s.dirs=[]);for(let t=0;t<n.length;t++){let[l,o,p,c=ks]=n[t];l&&(ds(l)&&(l={mounted:l,updated:l}),l.deep&&Fn(o),e.push({dir:l,instance:a,value:o,oldValue:void 0,arg:p,modifiers:c}))}return s}function ca(s,n,a,e){const t=s.dirs,l=n&&n.dirs;for(let o=0;o<t.length;o++){const p=t[o];l&&(p.oldValue=l[o].value);let c=p.dir[e];c&&(zn(),wn(c,a,8,[s.el,p,s,n]),qn())}}const bc=Symbol("_vte"),_c=s=>s.__isTeleport,Ga=s=>s&&(s.disabled||s.disabled===""),Jl=s=>s&&(s.defer||s.defer===""),Ql=s=>typeof SVGElement<"u"&&s instanceof SVGElement,Zl=s=>typeof MathMLElement=="function"&&s instanceof MathMLElement,Nt=(s,n)=>{const a=s&&s.to;return Ms(a)?n?n(a):null:a},yc={name:"Teleport",__isTeleport:!0,process(s,n,a,e,t,l,o,p,c,u){const{mc:r,pc:i,pbc:d,o:{insert:h,querySelector:j,createText:m,createComment:P}}=u,C=Ga(n.props);let{shapeFlag:b,children:x,dynamicChildren:f}=n;if(s==null){const k=n.el=m(""),w=n.anchor=m("");h(k,a,e),h(w,a,e);const y=(S,F)=>{b&16&&r(x,S,F,t,l,o,p,c)},M=()=>{const S=n.target=Nt(n.props,j),F=vc(S,n,m,h);S&&(o!=="svg"&&Ql(S)?o="svg":o!=="mathml"&&Zl(S)&&(o="mathml"),t&&t.isCE&&(t.ce._teleportTargets||(t.ce._teleportTargets=new Set)).add(S),C||(y(S,F),Ae(n,!1)))};C&&(y(a,w),Ae(n,!0)),Jl(n.props)?(n.el.__isMounted=!1,an(()=>{M(),delete n.el.__isMounted},l)):M()}else{if(Jl(n.props)&&s.el.__isMounted===!1){an(()=>{yc.process(s,n,a,e,t,l,o,p,c,u)},l);return}n.el=s.el,n.targetStart=s.targetStart;const k=n.anchor=s.anchor,w=n.target=s.target,y=n.targetAnchor=s.targetAnchor,M=Ga(s.props),S=M?a:w,F=M?k:y;if(o==="svg"||Ql(w)?o="svg":(o==="mathml"||Zl(w))&&(o="mathml"),f?(d(s.dynamicChildren,f,S,t,l,o,p),Al(s,n,!0)):c||i(s,n,S,F,t,l,o,p,!1),C)M?n.props&&s.props&&n.props.to!==s.props.to&&(n.props.to=s.props.to):ve(n,a,k,u,1);else if((n.props&&n.props.to)!==(s.props&&s.props.to)){const Z=n.target=Nt(n.props,j);Z&&ve(n,Z,null,u,0)}else M&&ve(n,w,y,u,1);Ae(n,C)}},remove(s,n,a,{um:e,o:{remove:t}},l){const{shapeFlag:o,children:p,anchor:c,targetStart:u,targetAnchor:r,target:i,props:d}=s;if(i&&(t(u),t(r)),l&&t(c),o&16){const h=l||!Ga(d);for(let j=0;j<p.length;j++){const m=p[j];e(m,n,a,h,!!m.dynamicChildren)}}},move:ve,hydrate:lu};function ve(s,n,a,{o:{insert:e},m:t},l=2){l===0&&e(s.targetAnchor,n,a);const{el:o,anchor:p,shapeFlag:c,children:u,props:r}=s,i=l===2;if(i&&e(o,n,a),(!i||Ga(r))&&c&16)for(let d=0;d<u.length;d++)t(u[d],n,a,2);i&&e(p,n,a)}function lu(s,n,a,e,t,l,{o:{nextSibling:o,parentNode:p,querySelector:c,insert:u,createText:r}},i){function d(m,P,C,b){P.anchor=i(o(m),P,p(m),a,e,t,l),P.targetStart=C,P.targetAnchor=b}const h=n.target=Nt(n.props,c),j=Ga(n.props);if(h){const m=h._lpa||h.firstChild;if(n.shapeFlag&16)if(j)d(s,n,m,m&&o(m));else{n.anchor=o(s);let P=m;for(;P;){if(P&&P.nodeType===8){if(P.data==="teleport start anchor")n.targetStart=P;else if(P.data==="teleport anchor"){n.targetAnchor=P,h._lpa=n.targetAnchor&&o(n.targetAnchor);break}}P=o(P)}n.targetAnchor||vc(h,n,r,u),i(m&&o(m),n,h,a,e,t,l)}Ae(n,j)}else j&&n.shapeFlag&16&&d(s,n,s,o(s));return n.anchor&&o(n.anchor)}const dv=yc;function Ae(s,n){const a=s.ctx;if(a&&a.ut){let e,t;for(n?(e=s.el,t=s.anchor):(e=s.targetStart,t=s.targetAnchor);e&&e!==t;)e.nodeType===1&&e.setAttribute("data-v-owner",a.uid),e=e.nextSibling;a.ut()}}function vc(s,n,a,e){const t=n.targetStart=a(""),l=n.targetAnchor=a("");return t[bc]=l,s&&(e(t,s),e(l,s)),l}const Ln=Symbol("_leaveCb"),we=Symbol("_enterCb");function ou(){const s={isMounted:!1,isLeaving:!1,isUnmounting:!1,leavingVNodes:new Map};return dn(()=>{s.isMounted=!0}),Tn(()=>{s.isUnmounting=!0}),s}const gn=[Function,Array],wc={mode:String,appear:Boolean,persisted:Boolean,onBeforeEnter:gn,onEnter:gn,onAfterEnter:gn,onEnterCancelled:gn,onBeforeLeave:gn,onLeave:gn,onAfterLeave:gn,onLeaveCancelled:gn,onBeforeAppear:gn,onAppear:gn,onAfterAppear:gn,onAppearCancelled:gn},xc=s=>{const n=s.subTree;return n.component?xc(n.component):n},pu={name:"BaseTransition",props:wc,setup(s,{slots:n}){const a=Gn(),e=ou();return()=>{const t=n.default&&Sc(n.default(),!0);if(!t||!t.length)return;const l=Cc(t),o=ys(s),{mode:p}=o;if(e.isLeaving)return bt(l);const c=so(l);if(!c)return bt(l);let u=Bt(c,o,e,a,i=>u=i);c.type!==ln&&ae(c,u);let r=a.subTree&&so(a.subTree);if(r&&r.type!==ln&&!da(r,c)&&xc(a).type!==ln){let i=Bt(r,o,e,a);if(ae(r,i),p==="out-in"&&c.type!==ln)return e.isLeaving=!0,i.afterLeave=()=>{e.isLeaving=!1,a.job.flags&8||a.update(),delete i.afterLeave,r=void 0},bt(l);p==="in-out"&&c.type!==ln?i.delayLeave=(d,h,j)=>{const m=kc(e,r);m[String(r.key)]=r,d[Ln]=()=>{h(),d[Ln]=void 0,delete u.delayedLeave,r=void 0},u.delayedLeave=()=>{j(),delete u.delayedLeave,r=void 0}}:r=void 0}else r&&(r=void 0);return l}}};function Cc(s){let n=s[0];if(s.length>1){for(const a of s)if(a.type!==ln){n=a;break}}return n}const cu=pu;function kc(s,n){const{leavingVNodes:a}=s;let e=a.get(n.type);return e||(e=Object.create(null),a.set(n.type,e)),e}function Bt(s,n,a,e,t){const{appear:l,mode:o,persisted:p=!1,onBeforeEnter:c,onEnter:u,onAfterEnter:r,onEnterCancelled:i,onBeforeLeave:d,onLeave:h,onAfterLeave:j,onLeaveCancelled:m,onBeforeAppear:P,onAppear:C,onAfterAppear:b,onAppearCancelled:x}=n,f=String(s.key),k=kc(a,s),w=(S,F)=>{S&&wn(S,e,9,F)},y=(S,F)=>{const Z=F[1];w(S,F),is(S)?S.every(B=>B.length<=1)&&Z():S.length<=1&&Z()},M={mode:o,persisted:p,beforeEnter(S){let F=c;if(!a.isMounted)if(l)F=P||c;else return;S[Ln]&&S[Ln](!0);const Z=k[f];Z&&da(s,Z)&&Z.el[Ln]&&Z.el[Ln](),w(F,[S])},enter(S){let F=u,Z=r,B=i;if(!a.isMounted)if(l)F=C||u,Z=b||r,B=x||i;else return;let O=!1;const V=S[we]=us=>{O||(O=!0,us?w(B,[S]):w(Z,[S]),M.delayedLeave&&M.delayedLeave(),S[we]=void 0)};F?y(F,[S,V]):V()},leave(S,F){const Z=String(s.key);if(S[we]&&S[we](!0),a.isUnmounting)return F();w(d,[S]);let B=!1;const O=S[Ln]=V=>{B||(B=!0,F(),V?w(m,[S]):w(j,[S]),S[Ln]=void 0,k[Z]===s&&delete k[Z])};k[Z]=s,h?y(h,[S,O]):O()},clone(S){const F=Bt(S,n,a,e,t);return t&&t(F),F}};return M}function bt(s){if(he(s))return s=aa(s),s.children=null,s}function so(s){if(!he(s))return _c(s.type)&&s.children?Cc(s.children):s;if(s.component)return s.component.subTree;const{shapeFlag:n,children:a}=s;if(a){if(n&16)return a[0];if(n&32&&ds(a.default))return a.default()}}function ae(s,n){s.shapeFlag&6&&s.component?(s.transition=n,ae(s.component.subTree,n)):s.shapeFlag&128?(s.ssContent.transition=n.clone(s.ssContent),s.ssFallback.transition=n.clone(s.ssFallback)):s.transition=n}function Sc(s,n=!1,a){let e=[],t=0;for(let l=0;l<s.length;l++){let o=s[l];const p=a==null?o.key:String(a)+String(o.key!=null?o.key:l);o.type===js?(o.patchFlag&128&&t++,e=e.concat(Sc(o.children,n,p))):(n||o.type!==ln)&&e.push(p!=null?aa(o,{key:p}):o)}if(t>1)for(let l=0;l<e.length;l++)e[l].patchFlag=-2;return e}function zs(s,n){return ds(s)?Gs({name:s.name},n,{setup:s}):s}function vl(s){s.ids=[s.ids[0]+s.ids[2]+++"-",0,0]}function Te(s){const n=Gn(),a=_l(null);if(n){const t=n.refs===ks?n.refs={}:n.refs;Object.defineProperty(t,s,{enumerable:!0,get:()=>a.value,set:l=>a.value=l})}return a}const ze=new WeakMap;function Ua(s,n,a,e,t=!1){if(is(s)){s.forEach((j,m)=>Ua(j,n&&(is(n)?n[m]:n),a,e,t));return}if(Wa(e)&&!t){e.shapeFlag&512&&e.type.__asyncResolved&&e.component.subTree.component&&Ua(s,n,a,e.component.subTree);return}const l=e.shapeFlag&4?tt(e.component):e.el,o=t?null:l,{i:p,r:c}=s,u=n&&n.r,r=p.refs===ks?p.refs={}:p.refs,i=p.setupState,d=ys(i),h=i===ks?Np:j=>xs(d,j);if(u!=null&&u!==c){if(no(n),Ms(u))r[u]=null,h(u)&&(i[u]=null);else if(Es(u)){u.value=null;const j=n;j.k&&(r[j.k]=null)}}if(ds(c))ie(c,p,12,[o,r]);else{const j=Ms(c),m=Es(c);if(j||m){const P=()=>{if(s.f){const C=j?h(c)?i[c]:r[c]:c.value;if(t)is(C)&&pl(C,l);else if(is(C))C.includes(l)||C.push(l);else if(j)r[c]=[l],h(c)&&(i[c]=r[c]);else{const b=[l];c.value=b,s.k&&(r[s.k]=b)}}else j?(r[c]=o,h(c)&&(i[c]=o)):m&&(c.value=o,s.k&&(r[s.k]=o))};if(o){const C=()=>{P(),ze.delete(s)};C.id=-1,ze.set(s,C),an(C,a)}else no(s),P()}}}function no(s){const n=ze.get(s);n&&(n.flags|=8,ze.delete(s))}const ao=s=>s.nodeType===8;st().requestIdleCallback;st().cancelIdleCallback;function ru(s,n){if(ao(s)&&s.data==="["){let a=1,e=s.nextSibling;for(;e;){if(e.nodeType===1){if(n(e)===!1)break}else if(ao(e))if(e.data==="]"){if(--a===0)break}else e.data==="["&&a++;e=e.nextSibling}}else n(s)}const Wa=s=>!!s.type.__asyncLoader;function iu(s){ds(s)&&(s={loader:s});const{loader:n,loadingComponent:a,errorComponent:e,delay:t=200,hydrate:l,timeout:o,suspensible:p=!0,onError:c}=s;let u=null,r,i=0;const d=()=>(i++,u=null,h()),h=()=>{let j;return u||(j=u=n().catch(m=>{if(m=m instanceof Error?m:new Error(String(m)),c)return new Promise((P,C)=>{c(m,()=>P(d()),()=>C(m),i+1)});throw m}).then(m=>j!==u&&u?u:(m&&(m.__esModule||m[Symbol.toStringTag]==="Module")&&(m=m.default),r=m,m)))};return zs({name:"AsyncComponentWrapper",__asyncLoader:h,__asyncHydrate(j,m,P){let C=!1;(m.bu||(m.bu=[])).push(()=>C=!0);const b=()=>{C||P()},x=l?()=>{const f=l(b,k=>ru(j,k));f&&(m.bum||(m.bum=[])).push(f)}:b;r?x():h().then(()=>!m.isUnmounted&&x())},get __asyncResolved(){return r},setup(){const j=Ys;if(vl(j),r)return()=>xe(r,j);const m=x=>{u=null,ue(x,j,13,!e)};if(p&&j.suspense||Ta)return h().then(x=>()=>xe(x,j)).catch(x=>(m(x),()=>e?bs(e,{error:x}):null));const P=ns(!1),C=ns(),b=ns(!!t);return t&&setTimeout(()=>{b.value=!1},t),o!=null&&setTimeout(()=>{if(!P.value&&!C.value){const x=new Error(`Async component timed out after ${o}ms.`);m(x),C.value=x}},o),h().then(()=>{P.value=!0,j.parent&&he(j.parent.vnode)&&j.parent.update()}).catch(x=>{m(x),C.value=x}),()=>{if(P.value&&r)return xe(r,j);if(C.value&&e)return bs(e,{error:C.value});if(a&&!b.value)return xe(a,j)}}})}function xe(s,n){const{ref:a,props:e,children:t,ce:l}=n.vnode,o=bs(s,e,t);return o.ref=a,o.ce=l,delete n.vnode.ce,o}const he=s=>s.type.__isKeepAlive;function Pc(s,n){Tc(s,"a",n)}function Ac(s,n){Tc(s,"da",n)}function Tc(s,n,a=Ys){const e=s.__wdc||(s.__wdc=()=>{let t=a;for(;t;){if(t.isDeactivated)return;t=t.parent}return s()});if(at(n,e,a),a){let t=a.parent;for(;t&&t.parent;)he(t.parent.vnode)&&uu(e,n,a,t),t=t.parent}}function uu(s,n,a,e){const t=at(n,s,e,!0);wl(()=>{pl(e[n],t)},a)}function at(s,n,a=Ys,e=!1){if(a){const t=a[s]||(a[s]=[]),l=n.__weh||(n.__weh=(...o)=>{zn();const p=je(a),c=wn(n,a,s,o);return p(),qn(),c});return e?t.unshift(l):t.push(l),l}}const Vn=s=>(n,a=Ys)=>{(!Ta||s==="sp")&&at(s,(...e)=>n(...e),a)},hu=Vn("bm"),dn=Vn("m"),du=Vn("bu"),mu=Vn("u"),Tn=Vn("bum"),wl=Vn("um"),ju=Vn("sp"),gu=Vn("rtg"),fu=Vn("rtc");function bu(s,n=Ys){at("ec",s,n)}const xl="components",_u="directives";function de(s,n){return kl(xl,s,!0,n)||s}const Rc=Symbol.for("v-ndc");function yu(s){return Ms(s)?kl(xl,s,!1)||s:s||Rc}function Cl(s){return kl(_u,s)}function kl(s,n,a=!0,e=!1){const t=cn||Ys;if(t){const l=t.type;if(s===xl){const p=ch(l,!1);if(p&&(p===n||p===yn(n)||p===Ze(yn(n))))return l}const o=eo(t[s]||l[s],n)||eo(t.appContext[s],n);return!o&&e?l:o}}function eo(s,n){return s&&(s[n]||s[yn(n)]||s[Ze(yn(n))])}function Ts(s,n,a,e){let t;const l=a,o=is(s);if(o||Ms(s)){const p=o&&sa(s);let c=!1,u=!1;p&&(c=!bn(s),u=na(s),s=nt(s)),t=new Array(s.length);for(let r=0,i=s.length;r<i;r++)t[r]=n(c?u?$e(Ws(s[r])):Ws(s[r]):s[r],r,void 0,l)}else if(typeof s=="number"){t=new Array(s);for(let p=0;p<s;p++)t[p]=n(p+1,p,void 0,l)}else if(Rs(s))if(s[Symbol.iterator])t=Array.from(s,(p,c)=>n(p,c,void 0,l));else{const p=Object.keys(s);t=new Array(p.length);for(let c=0,u=p.length;c<u;c++){const r=p[c];t[c]=n(s[r],r,c,l)}}else t=[];return t}const zt=s=>s?Yc(s)?tt(s):zt(s.parent):null,Ka=Gs(Object.create(null),{$:s=>s,$el:s=>s.vnode.el,$data:s=>s.data,$props:s=>s.props,$attrs:s=>s.attrs,$slots:s=>s.slots,$refs:s=>s.refs,$parent:s=>zt(s.parent),$root:s=>zt(s.root),$host:s=>s.ce,$emit:s=>s.emit,$options:s=>Mc(s),$forceUpdate:s=>s.f||(s.f=()=>{yl(s.update)}),$nextTick:s=>s.n||(s.n=Ks.bind(s.proxy)),$watch:s=>qu.bind(s)}),_t=(s,n)=>s!==ks&&!s.__isScriptSetup&&xs(s,n),vu={get({_:s},n){if(n==="__v_skip")return!0;const{ctx:a,setupState:e,data:t,props:l,accessCache:o,type:p,appContext:c}=s;let u;if(n[0]!=="$"){const h=o[n];if(h!==void 0)switch(h){case 1:return e[n];case 2:return t[n];case 4:return a[n];case 3:return l[n]}else{if(_t(e,n))return o[n]=1,e[n];if(t!==ks&&xs(t,n))return o[n]=2,t[n];if((u=s.propsOptions[0])&&xs(u,n))return o[n]=3,l[n];if(a!==ks&&xs(a,n))return o[n]=4,a[n];qt&&(o[n]=0)}}const r=Ka[n];let i,d;if(r)return n==="$attrs"&&Qs(s.attrs,"get",""),r(s);if((i=p.__cssModules)&&(i=i[n]))return i;if(a!==ks&&xs(a,n))return o[n]=4,a[n];if(d=c.config.globalProperties,xs(d,n))return d[n]},set({_:s},n,a){const{data:e,setupState:t,ctx:l}=s;return _t(t,n)?(t[n]=a,!0):e!==ks&&xs(e,n)?(e[n]=a,!0):xs(s.props,n)||n[0]==="$"&&n.slice(1)in s?!1:(l[n]=a,!0)},has({_:{data:s,setupState:n,accessCache:a,ctx:e,appContext:t,propsOptions:l,type:o}},p){let c,u;return!!(a[p]||s!==ks&&p[0]!=="$"&&xs(s,p)||_t(n,p)||(c=l[0])&&xs(c,p)||xs(e,p)||xs(Ka,p)||xs(t.config.globalProperties,p)||(u=o.__cssModules)&&u[p])},defineProperty(s,n,a){return a.get!=null?s._.accessCache[n]=0:xs(a,"value")&&this.set(s,n,a.value,null),Reflect.defineProperty(s,n,a)}};function to(s){return is(s)?s.reduce((n,a)=>(n[a]=null,n),{}):s}let qt=!0;function wu(s){const n=Mc(s),a=s.proxy,e=s.ctx;qt=!1,n.beforeCreate&&lo(n.beforeCreate,s,"bc");const{data:t,computed:l,methods:o,watch:p,provide:c,inject:u,created:r,beforeMount:i,mounted:d,beforeUpdate:h,updated:j,activated:m,deactivated:P,beforeDestroy:C,beforeUnmount:b,destroyed:x,unmounted:f,render:k,renderTracked:w,renderTriggered:y,errorCaptured:M,serverPrefetch:S,expose:F,inheritAttrs:Z,components:B,directives:O,filters:V}=n;if(u&&xu(u,e,null),o)for(const Y in o){const hs=o[Y];ds(hs)&&(e[Y]=hs.bind(a))}if(t){const Y=t.call(a,a);Rs(Y)&&(s.data=re(Y))}if(qt=!0,l)for(const Y in l){const hs=l[Y],Fs=ds(hs)?hs.bind(a,a):ds(hs.get)?hs.get.bind(a,a):An,$s=!ds(hs)&&ds(hs.set)?hs.set.bind(a):An,Hs=ss({get:Fs,set:$s});Object.defineProperty(e,Y,{enumerable:!0,configurable:!0,get:()=>Hs.value,set:Vs=>Hs.value=Vs})}if(p)for(const Y in p)Ec(p[Y],e,a,Y);if(c){const Y=ds(c)?c.call(a):c;Reflect.ownKeys(Y).forEach(hs=>{Re(hs,Y[hs])})}r&&lo(r,s,"c");function os(Y,hs){is(hs)?hs.forEach(Fs=>Y(Fs.bind(a))):hs&&Y(hs.bind(a))}if(os(hu,i),os(dn,d),os(du,h),os(mu,j),os(Pc,m),os(Ac,P),os(bu,M),os(fu,w),os(gu,y),os(Tn,b),os(wl,f),os(ju,S),is(F))if(F.length){const Y=s.exposed||(s.exposed={});F.forEach(hs=>{Object.defineProperty(Y,hs,{get:()=>a[hs],set:Fs=>a[hs]=Fs,enumerable:!0})})}else s.exposed||(s.exposed={});k&&s.render===An&&(s.render=k),Z!=null&&(s.inheritAttrs=Z),B&&(s.components=B),O&&(s.directives=O),S&&vl(s)}function xu(s,n,a=An){is(s)&&(s=Ht(s));for(const e in s){const t=s[e];let l;Rs(t)?"default"in t?l=on(t.from||e,t.default,!0):l=on(t.from||e):l=on(t),Es(l)?Object.defineProperty(n,e,{enumerable:!0,configurable:!0,get:()=>l.value,set:o=>l.value=o}):n[e]=l}}function lo(s,n,a){wn(is(s)?s.map(e=>e.bind(n.proxy)):s.bind(n.proxy),n,a)}function Ec(s,n,a,e){let t=e.includes(".")?Vc(a,e):()=>a[e];if(Ms(s)){const l=n[s];ds(l)&&Bs(t,l)}else if(ds(s))Bs(t,s.bind(a));else if(Rs(s))if(is(s))s.forEach(l=>Ec(l,n,a,e));else{const l=ds(s.handler)?s.handler.bind(a):n[s.handler];ds(l)&&Bs(t,l,s)}}function Mc(s){const n=s.type,{mixins:a,extends:e}=n,{mixins:t,optionsCache:l,config:{optionMergeStrategies:o}}=s.appContext,p=l.get(n);let c;return p?c=p:!t.length&&!a&&!e?c=n:(c={},t.length&&t.forEach(u=>qe(c,u,o,!0)),qe(c,n,o)),Rs(n)&&l.set(n,c),c}function qe(s,n,a,e=!1){const{mixins:t,extends:l}=n;l&&qe(s,l,a,!0),t&&t.forEach(o=>qe(s,o,a,!0));for(const o in n)if(!(e&&o==="expose")){const p=Cu[o]||a&&a[o];s[o]=p?p(s[o],n[o]):n[o]}return s}const Cu={data:oo,props:po,emits:po,methods:za,computed:za,beforeCreate:nn,created:nn,beforeMount:nn,mounted:nn,beforeUpdate:nn,updated:nn,beforeDestroy:nn,beforeUnmount:nn,destroyed:nn,unmounted:nn,activated:nn,deactivated:nn,errorCaptured:nn,serverPrefetch:nn,components:za,directives:za,watch:Su,provide:oo,inject:ku};function oo(s,n){return n?s?function(){return Gs(ds(s)?s.call(this,this):s,ds(n)?n.call(this,this):n)}:n:s}function ku(s,n){return za(Ht(s),Ht(n))}function Ht(s){if(is(s)){const n={};for(let a=0;a<s.length;a++)n[s[a]]=s[a];return n}return s}function nn(s,n){return s?[...new Set([].concat(s,n))]:n}function za(s,n){return s?Gs(Object.create(null),s,n):n}function po(s,n){return s?is(s)&&is(n)?[...new Set([...s,...n])]:Gs(Object.create(null),to(s),to(n??{})):n}function Su(s,n){if(!s)return n;if(!n)return s;const a=Gs(Object.create(null),s);for(const e in n)a[e]=nn(s[e],n[e]);return a}function Dc(){return{app:null,config:{isNativeTag:Np,performance:!1,globalProperties:{},optionMergeStrategies:{},errorHandler:void 0,warnHandler:void 0,compilerOptions:{}},mixins:[],components:{},directives:{},provides:Object.create(null),optionsCache:new WeakMap,propsCache:new WeakMap,emitsCache:new WeakMap}}let Pu=0;function Au(s,n){return function(e,t=null){ds(e)||(e=Gs({},e)),t!=null&&!Rs(t)&&(t=null);const l=Dc(),o=new WeakSet,p=[];let c=!1;const u=l.app={_uid:Pu++,_component:e,_props:t,_container:null,_context:l,_instance:null,version:ih,get config(){return l.config},set config(r){},use(r,...i){return o.has(r)||(r&&ds(r.install)?(o.add(r),r.install(u,...i)):ds(r)&&(o.add(r),r(u,...i))),u},mixin(r){return l.mixins.includes(r)||l.mixins.push(r),u},component(r,i){return i?(l.components[r]=i,u):l.components[r]},directive(r,i){return i?(l.directives[r]=i,u):l.directives[r]},mount(r,i,d){if(!c){const h=u._ceVNode||bs(e,t);return h.appContext=l,d===!0?d="svg":d===!1&&(d=void 0),s(h,r,d),c=!0,u._container=r,r.__vue_app__=u,tt(h.component)}},onUnmount(r){p.push(r)},unmount(){c&&(wn(p,u._instance,16),s(null,u._container),delete u._container.__vue_app__)},provide(r,i){return l.provides[r]=i,u},runWithContext(r){const i=ga;ga=u;try{return r()}finally{ga=i}}};return u}}let ga=null;function Re(s,n){if(Ys){let a=Ys.provides;const e=Ys.parent&&Ys.parent.provides;e===a&&(a=Ys.provides=Object.create(e)),a[s]=n}}function on(s,n,a=!1){const e=Gn();if(e||ga){let t=ga?ga._context.provides:e?e.parent==null||e.ce?e.vnode.appContext&&e.vnode.appContext.provides:e.parent.provides:void 0;if(t&&s in t)return t[s];if(arguments.length>1)return a&&ds(n)?n.call(e&&e.proxy):n}}function Lc(){return!!(Gn()||ga)}const Oc={},Fc=()=>Object.create(Oc),$c=s=>Object.getPrototypeOf(s)===Oc;function Tu(s,n,a,e=!1){const t={},l=Fc();s.propsDefaults=Object.create(null),Ic(s,n,t,l);for(const o in s.propsOptions[0])o in t||(t[o]=void 0);a?s.props=e?t:ic(t):s.type.props?s.props=t:s.props=l,s.attrs=l}function Ru(s,n,a,e){const{props:t,attrs:l,vnode:{patchFlag:o}}=s,p=ys(t),[c]=s.propsOptions;let u=!1;if((e||o>0)&&!(o&16)){if(o&8){const r=s.vnode.dynamicProps;for(let i=0;i<r.length;i++){let d=r[i];if(et(s.emitsOptions,d))continue;const h=n[d];if(c)if(xs(l,d))h!==l[d]&&(l[d]=h,u=!0);else{const j=yn(d);t[j]=Vt(c,p,j,h,s,!1)}else h!==l[d]&&(l[d]=h,u=!0)}}}else{Ic(s,n,t,l)&&(u=!0);let r;for(const i in p)(!n||!xs(n,i)&&((r=la(i))===i||!xs(n,r)))&&(c?a&&(a[i]!==void 0||a[r]!==void 0)&&(t[i]=Vt(c,p,i,void 0,s,!0)):delete t[i]);if(l!==p)for(const i in l)(!n||!xs(n,i))&&(delete l[i],u=!0)}u&&On(s.attrs,"set","")}function Ic(s,n,a,e){const[t,l]=s.propsOptions;let o=!1,p;if(n)for(let c in n){if(qa(c))continue;const u=n[c];let r;t&&xs(t,r=yn(c))?!l||!l.includes(r)?a[r]=u:(p||(p={}))[r]=u:et(s.emitsOptions,c)||(!(c in e)||u!==e[c])&&(e[c]=u,o=!0)}if(l){const c=ys(a),u=p||ks;for(let r=0;r<l.length;r++){const i=l[r];a[i]=Vt(t,c,i,u[i],s,!xs(u,i))}}return o}function Vt(s,n,a,e,t,l){const o=s[a];if(o!=null){const p=xs(o,"default");if(p&&e===void 0){const c=o.default;if(o.type!==Function&&!o.skipFactory&&ds(c)){const{propsDefaults:u}=t;if(a in u)e=u[a];else{const r=je(t);e=u[a]=c.call(null,n),r()}}else e=c;t.ce&&t.ce._setProp(a,e)}o[0]&&(l&&!p?e=!1:o[1]&&(e===""||e===la(a))&&(e=!0))}return e}const Eu=new WeakMap;function Nc(s,n,a=!1){const e=a?Eu:n.propsCache,t=e.get(s);if(t)return t;const l=s.props,o={},p=[];let c=!1;if(!ds(s)){const r=i=>{c=!0;const[d,h]=Nc(i,n,!0);Gs(o,d),h&&p.push(...h)};!a&&n.mixins.length&&n.mixins.forEach(r),s.extends&&r(s.extends),s.mixins&&s.mixins.forEach(r)}if(!l&&!c)return Rs(s)&&e.set(s,Ca),Ca;if(is(l))for(let r=0;r<l.length;r++){const i=yn(l[r]);co(i)&&(o[i]=ks)}else if(l)for(const r in l){const i=yn(r);if(co(i)){const d=l[r],h=o[i]=is(d)||ds(d)?{type:d}:Gs({},d),j=h.type;let m=!1,P=!0;if(is(j))for(let C=0;C<j.length;++C){const b=j[C],x=ds(b)&&b.name;if(x==="Boolean"){m=!0;break}else x==="String"&&(P=!1)}else m=ds(j)&&j.name==="Boolean";h[0]=m,h[1]=P,(m||xs(h,"default"))&&p.push(i)}}const u=[o,p];return Rs(s)&&e.set(s,u),u}function co(s){return s[0]!=="$"&&!qa(s)}const Sl=s=>s==="_"||s==="_ctx"||s==="$stable",Pl=s=>is(s)?s.map(Pn):[Pn(s)],Mu=(s,n,a)=>{if(n._n)return n;const e=hn((...t)=>Pl(n(...t)),a);return e._c=!1,e},Bc=(s,n,a)=>{const e=s._ctx;for(const t in s){if(Sl(t))continue;const l=s[t];if(ds(l))n[t]=Mu(t,l,e);else if(l!=null){const o=Pl(l);n[t]=()=>o}}},zc=(s,n)=>{const a=Pl(n);s.slots.default=()=>a},qc=(s,n,a)=>{for(const e in n)(a||!Sl(e))&&(s[e]=n[e])},Du=(s,n,a)=>{const e=s.slots=Fc();if(s.vnode.shapeFlag&32){const t=n._;t?(qc(e,n,a),a&&Vp(e,"_",t,!0)):Bc(n,e)}else n&&zc(s,n)},Lu=(s,n,a)=>{const{vnode:e,slots:t}=s;let l=!0,o=ks;if(e.shapeFlag&32){const p=n._;p?a&&p===1?l=!1:qc(t,n,a):(l=!n.$stable,Bc(n,t)),o=n}else n&&(zc(s,n),o={default:1});if(l)for(const p in t)!Sl(p)&&o[p]==null&&delete t[p]},an=Xu;function Ou(s){return Fu(s)}function Fu(s,n){const a=st();a.__VUE__=!0;const{insert:e,remove:t,patchProp:l,createElement:o,createText:p,createComment:c,setText:u,setElementText:r,parentNode:i,nextSibling:d,setScopeId:h=An,insertStaticContent:j}=s,m=(g,_,A,$=null,N=null,I=null,K=void 0,W=null,U=!!_.dynamicChildren)=>{if(g===_)return;g&&!da(g,_)&&($=R(g),Vs(g,N,I,!0),g=null),_.patchFlag===-2&&(U=!1,_.dynamicChildren=null);const{type:z,ref:ps,shapeFlag:J}=_;switch(z){case me:P(g,_,A,$);break;case ln:C(g,_,A,$);break;case Ee:g==null&&b(_,A,$,K);break;case js:B(g,_,A,$,N,I,K,W,U);break;default:J&1?k(g,_,A,$,N,I,K,W,U):J&6?O(g,_,A,$,N,I,K,W,U):(J&64||J&128)&&z.process(g,_,A,$,N,I,K,W,U,X)}ps!=null&&N?Ua(ps,g&&g.ref,I,_||g,!_):ps==null&&g&&g.ref!=null&&Ua(g.ref,null,I,g,!0)},P=(g,_,A,$)=>{if(g==null)e(_.el=p(_.children),A,$);else{const N=_.el=g.el;_.children!==g.children&&u(N,_.children)}},C=(g,_,A,$)=>{g==null?e(_.el=c(_.children||""),A,$):_.el=g.el},b=(g,_,A,$)=>{[g.el,g.anchor]=j(g.children,_,A,$,g.el,g.anchor)},x=({el:g,anchor:_},A,$)=>{let N;for(;g&&g!==_;)N=d(g),e(g,A,$),g=N;e(_,A,$)},f=({el:g,anchor:_})=>{let A;for(;g&&g!==_;)A=d(g),t(g),g=A;t(_)},k=(g,_,A,$,N,I,K,W,U)=>{if(_.type==="svg"?K="svg":_.type==="math"&&(K="mathml"),g==null)w(_,A,$,N,I,K,W,U);else{const z=g.el&&g.el._isVueCE?g.el:null;try{z&&z._beginPatch(),S(g,_,N,I,K,W,U)}finally{z&&z._endPatch()}}},w=(g,_,A,$,N,I,K,W)=>{let U,z;const{props:ps,shapeFlag:J,transition:Q,dirs:T}=g;if(U=g.el=o(g.type,I,ps&&ps.is,ps),J&8?r(U,g.children):J&16&&M(g.children,U,null,$,N,yt(g,I),K,W),T&&ca(g,null,$,"created"),y(U,g,g.scopeId,K,$),ps){for(const ts in ps)ts!=="value"&&!qa(ts)&&l(U,ts,null,ps[ts],I,$);"value"in ps&&l(U,"value",null,ps.value,I),(z=ps.onVnodeBeforeMount)&&kn(z,$,g)}T&&ca(g,null,$,"beforeMount");const D=$u(N,Q);D&&Q.beforeEnter(U),e(U,_,A),((z=ps&&ps.onVnodeMounted)||D||T)&&an(()=>{z&&kn(z,$,g),D&&Q.enter(U),T&&ca(g,null,$,"mounted")},N)},y=(g,_,A,$,N)=>{if(A&&h(g,A),$)for(let I=0;I<$.length;I++)h(g,$[I]);if(N){let I=N.subTree;if(_===I||Uc(I.type)&&(I.ssContent===_||I.ssFallback===_)){const K=N.vnode;y(g,K,K.scopeId,K.slotScopeIds,N.parent)}}},M=(g,_,A,$,N,I,K,W,U=0)=>{for(let z=U;z<g.length;z++){const ps=g[z]=W?Xn(g[z]):Pn(g[z]);m(null,ps,_,A,$,N,I,K,W)}},S=(g,_,A,$,N,I,K)=>{const W=_.el=g.el;let{patchFlag:U,dynamicChildren:z,dirs:ps}=_;U|=g.patchFlag&16;const J=g.props||ks,Q=_.props||ks;let T;if(A&&ra(A,!1),(T=Q.onVnodeBeforeUpdate)&&kn(T,A,_,g),ps&&ca(_,g,A,"beforeUpdate"),A&&ra(A,!0),(J.innerHTML&&Q.innerHTML==null||J.textContent&&Q.textContent==null)&&r(W,""),z?F(g.dynamicChildren,z,W,A,$,yt(_,N),I):K||hs(g,_,W,null,A,$,yt(_,N),I,!1),U>0){if(U&16)Z(W,J,Q,A,N);else if(U&2&&J.class!==Q.class&&l(W,"class",null,Q.class,N),U&4&&l(W,"style",J.style,Q.style,N),U&8){const D=_.dynamicProps;for(let ts=0;ts<D.length;ts++){const cs=D[ts],Ds=J[cs],qs=Q[cs];(qs!==Ds||cs==="value")&&l(W,cs,Ds,qs,N,A)}}U&1&&g.children!==_.children&&r(W,_.children)}else!K&&z==null&&Z(W,J,Q,A,N);((T=Q.onVnodeUpdated)||ps)&&an(()=>{T&&kn(T,A,_,g),ps&&ca(_,g,A,"updated")},$)},F=(g,_,A,$,N,I,K)=>{for(let W=0;W<_.length;W++){const U=g[W],z=_[W],ps=U.el&&(U.type===js||!da(U,z)||U.shapeFlag&198)?i(U.el):A;m(U,z,ps,null,$,N,I,K,!0)}},Z=(g,_,A,$,N)=>{if(_!==A){if(_!==ks)for(const I in _)!qa(I)&&!(I in A)&&l(g,I,_[I],null,N,$);for(const I in A){if(qa(I))continue;const K=A[I],W=_[I];K!==W&&I!=="value"&&l(g,I,W,K,N,$)}"value"in A&&l(g,"value",_.value,A.value,N)}},B=(g,_,A,$,N,I,K,W,U)=>{const z=_.el=g?g.el:p(""),ps=_.anchor=g?g.anchor:p("");let{patchFlag:J,dynamicChildren:Q,slotScopeIds:T}=_;T&&(W=W?W.concat(T):T),g==null?(e(z,A,$),e(ps,A,$),M(_.children||[],A,ps,N,I,K,W,U)):J>0&&J&64&&Q&&g.dynamicChildren?(F(g.dynamicChildren,Q,A,N,I,K,W),(_.key!=null||N&&_===N.subTree)&&Al(g,_,!0)):hs(g,_,A,ps,N,I,K,W,U)},O=(g,_,A,$,N,I,K,W,U)=>{_.slotScopeIds=W,g==null?_.shapeFlag&512?N.ctx.activate(_,A,$,K,U):V(_,A,$,N,I,K,U):us(g,_,U)},V=(g,_,A,$,N,I,K)=>{const W=g.component=eh(g,$,N);if(he(g)&&(W.ctx.renderer=X),th(W,!1,K),W.asyncDep){if(N&&N.registerDep(W,os,K),!g.el){const U=W.subTree=bs(ln);C(null,U,_,A),g.placeholder=U.el}}else os(W,g,_,A,N,I,K)},us=(g,_,A)=>{const $=_.component=g.component;if(Ku(g,_,A))if($.asyncDep&&!$.asyncResolved){Y($,_,A);return}else $.next=_,$.update();else _.el=g.el,$.vnode=_},os=(g,_,A,$,N,I,K)=>{const W=()=>{if(g.isMounted){let{next:J,bu:Q,u:T,parent:D,vnode:ts}=g;{const jn=Hc(g);if(jn){J&&(J.el=ts.el,Y(g,J,K)),jn.asyncDep.then(()=>{g.isUnmounted||W()});return}}let cs=J,Ds;ra(g,!1),J?(J.el=ts.el,Y(g,J,K)):J=ts,Q&&Pe(Q),(Ds=J.props&&J.props.onVnodeBeforeUpdate)&&kn(Ds,D,J,ts),ra(g,!0);const qs=io(g),Us=g.subTree;g.subTree=qs,m(Us,qs,i(Us.el),R(Us),g,N,I),J.el=qs.el,cs===null&&Yu(g,qs.el),T&&an(T,N),(Ds=J.props&&J.props.onVnodeUpdated)&&an(()=>kn(Ds,D,J,ts),N)}else{let J;const{el:Q,props:T}=_,{bm:D,m:ts,parent:cs,root:Ds,type:qs}=g,Us=Wa(_);ra(g,!1),D&&Pe(D),!Us&&(J=T&&T.onVnodeBeforeMount)&&kn(J,cs,_),ra(g,!0);{Ds.ce&&Ds.ce._def.shadowRoot!==!1&&Ds.ce._injectChildStyle(qs);const jn=g.subTree=io(g);m(null,jn,A,$,g,N,I),_.el=jn.el}if(ts&&an(ts,N),!Us&&(J=T&&T.onVnodeMounted)){const jn=_;an(()=>kn(J,cs,jn),N)}(_.shapeFlag&256||cs&&Wa(cs.vnode)&&cs.vnode.shapeFlag&256)&&g.a&&an(g.a,N),g.isMounted=!0,_=A=$=null}};g.scope.on();const U=g.effect=new Xp(W);g.scope.off();const z=g.update=U.run.bind(U),ps=g.job=U.runIfDirty.bind(U);ps.i=g,ps.id=g.uid,U.scheduler=()=>yl(ps),ra(g,!0),z()},Y=(g,_,A)=>{_.component=g;const $=g.vnode.props;g.vnode=_,g.next=null,Ru(g,_.props,$,A),Lu(g,_.children,A),zn(),Xl(g),qn()},hs=(g,_,A,$,N,I,K,W,U=!1)=>{const z=g&&g.children,ps=g?g.shapeFlag:0,J=_.children,{patchFlag:Q,shapeFlag:T}=_;if(Q>0){if(Q&128){$s(z,J,A,$,N,I,K,W,U);return}else if(Q&256){Fs(z,J,A,$,N,I,K,W,U);return}}T&8?(ps&16&&fs(z,N,I),J!==z&&r(A,J)):ps&16?T&16?$s(z,J,A,$,N,I,K,W,U):fs(z,N,I,!0):(ps&8&&r(A,""),T&16&&M(J,A,$,N,I,K,W,U))},Fs=(g,_,A,$,N,I,K,W,U)=>{g=g||Ca,_=_||Ca;const z=g.length,ps=_.length,J=Math.min(z,ps);let Q;for(Q=0;Q<J;Q++){const T=_[Q]=U?Xn(_[Q]):Pn(_[Q]);m(g[Q],T,A,null,N,I,K,W,U)}z>ps?fs(g,N,I,!0,!1,J):M(_,A,$,N,I,K,W,U,J)},$s=(g,_,A,$,N,I,K,W,U)=>{let z=0;const ps=_.length;let J=g.length-1,Q=ps-1;for(;z<=J&&z<=Q;){const T=g[z],D=_[z]=U?Xn(_[z]):Pn(_[z]);if(da(T,D))m(T,D,A,null,N,I,K,W,U);else break;z++}for(;z<=J&&z<=Q;){const T=g[J],D=_[Q]=U?Xn(_[Q]):Pn(_[Q]);if(da(T,D))m(T,D,A,null,N,I,K,W,U);else break;J--,Q--}if(z>J){if(z<=Q){const T=Q+1,D=T<ps?_[T].el:$;for(;z<=Q;)m(null,_[z]=U?Xn(_[z]):Pn(_[z]),A,D,N,I,K,W,U),z++}}else if(z>Q)for(;z<=J;)Vs(g[z],N,I,!0),z++;else{const T=z,D=z,ts=new Map;for(z=D;z<=Q;z++){const pn=_[z]=U?Xn(_[z]):Pn(_[z]);pn.key!=null&&ts.set(pn.key,z)}let cs,Ds=0;const qs=Q-D+1;let Us=!1,jn=0;const ba=new Array(qs);for(z=0;z<qs;z++)ba[z]=0;for(z=T;z<=J;z++){const pn=g[z];if(Ds>=qs){Vs(pn,N,I,!0);continue}let Cn;if(pn.key!=null)Cn=ts.get(pn.key);else for(cs=D;cs<=Q;cs++)if(ba[cs-D]===0&&da(pn,_[cs])){Cn=cs;break}Cn===void 0?Vs(pn,N,I,!0):(ba[Cn-D]=z+1,Cn>=jn?jn=Cn:Us=!0,m(pn,_[Cn],A,null,N,I,K,W,U),Ds++)}const Nl=Us?Iu(ba):Ca;for(cs=Nl.length-1,z=qs-1;z>=0;z--){const pn=D+z,Cn=_[pn],Bl=_[pn+1],zl=pn+1<ps?Bl.el||Bl.placeholder:$;ba[z]===0?m(null,Cn,A,zl,N,I,K,W,U):Us&&(cs<0||z!==Nl[cs]?Hs(Cn,A,zl,2):cs--)}}},Hs=(g,_,A,$,N=null)=>{const{el:I,type:K,transition:W,children:U,shapeFlag:z}=g;if(z&6){Hs(g.component.subTree,_,A,$);return}if(z&128){g.suspense.move(_,A,$);return}if(z&64){K.move(g,_,A,X);return}if(K===js){e(I,_,A);for(let J=0;J<U.length;J++)Hs(U[J],_,A,$);e(g.anchor,_,A);return}if(K===Ee){x(g,_,A);return}if($!==2&&z&1&&W)if($===0)W.beforeEnter(I),e(I,_,A),an(()=>W.enter(I),N);else{const{leave:J,delayLeave:Q,afterLeave:T}=W,D=()=>{g.ctx.isUnmounted?t(I):e(I,_,A)},ts=()=>{I._isLeaving&&I[Ln](!0),J(I,()=>{D(),T&&T()})};Q?Q(I,D,ts):ts()}else e(I,_,A)},Vs=(g,_,A,$=!1,N=!1)=>{const{type:I,props:K,ref:W,children:U,dynamicChildren:z,shapeFlag:ps,patchFlag:J,dirs:Q,cacheIndex:T}=g;if(J===-2&&(N=!1),W!=null&&(zn(),Ua(W,null,A,g,!0),qn()),T!=null&&(_.renderCache[T]=void 0),ps&256){_.ctx.deactivate(g);return}const D=ps&1&&Q,ts=!Wa(g);let cs;if(ts&&(cs=K&&K.onVnodeBeforeUnmount)&&kn(cs,_,g),ps&6)_s(g.component,A,$);else{if(ps&128){g.suspense.unmount(A,$);return}D&&ca(g,null,_,"beforeUnmount"),ps&64?g.type.remove(g,_,A,X,$):z&&!z.hasOnce&&(I!==js||J>0&&J&64)?fs(z,_,A,!1,!0):(I===js&&J&384||!N&&ps&16)&&fs(U,_,A),$&&es(g)}(ts&&(cs=K&&K.onVnodeUnmounted)||D)&&an(()=>{cs&&kn(cs,_,g),D&&ca(g,null,_,"unmounted")},A)},es=g=>{const{type:_,el:A,anchor:$,transition:N}=g;if(_===js){ms(A,$);return}if(_===Ee){f(g);return}const I=()=>{t(A),N&&!N.persisted&&N.afterLeave&&N.afterLeave()};if(g.shapeFlag&1&&N&&!N.persisted){const{leave:K,delayLeave:W}=N,U=()=>K(A,I);W?W(g.el,I,U):U()}else I()},ms=(g,_)=>{let A;for(;g!==_;)A=d(g),t(g),g=A;t(_)},_s=(g,_,A)=>{const{bum:$,scope:N,job:I,subTree:K,um:W,m:U,a:z}=g;ro(U),ro(z),$&&Pe($),N.stop(),I&&(I.flags|=8,Vs(K,g,_,A)),W&&an(W,_),an(()=>{g.isUnmounted=!0},_)},fs=(g,_,A,$=!1,N=!1,I=0)=>{for(let K=I;K<g.length;K++)Vs(g[K],_,A,$,N)},R=g=>{if(g.shapeFlag&6)return R(g.component.subTree);if(g.shapeFlag&128)return g.suspense.next();const _=d(g.anchor||g.el),A=_&&_[bc];return A?d(A):_};let H=!1;const q=(g,_,A)=>{g==null?_._vnode&&Vs(_._vnode,null,null,!0):m(_._vnode||null,g,_,null,null,null,A),_._vnode=g,H||(H=!0,Xl(),jc(),H=!1)},X={p:m,um:Vs,m:Hs,r:es,mt:V,mc:M,pc:hs,pbc:F,n:R,o:s};return{render:q,hydrate:void 0,createApp:Au(q)}}function yt({type:s,props:n},a){return a==="svg"&&s==="foreignObject"||a==="mathml"&&s==="annotation-xml"&&n&&n.encoding&&n.encoding.includes("html")?void 0:a}function ra({effect:s,job:n},a){a?(s.flags|=32,n.flags|=4):(s.flags&=-33,n.flags&=-5)}function $u(s,n){return(!s||s&&!s.pendingBranch)&&n&&!n.persisted}function Al(s,n,a=!1){const e=s.children,t=n.children;if(is(e)&&is(t))for(let l=0;l<e.length;l++){const o=e[l];let p=t[l];p.shapeFlag&1&&!p.dynamicChildren&&((p.patchFlag<=0||p.patchFlag===32)&&(p=t[l]=Xn(t[l]),p.el=o.el),!a&&p.patchFlag!==-2&&Al(o,p)),p.type===me&&p.patchFlag!==-1&&(p.el=o.el),p.type===ln&&!p.el&&(p.el=o.el)}}function Iu(s){const n=s.slice(),a=[0];let e,t,l,o,p;const c=s.length;for(e=0;e<c;e++){const u=s[e];if(u!==0){if(t=a[a.length-1],s[t]<u){n[e]=t,a.push(e);continue}for(l=0,o=a.length-1;l<o;)p=l+o>>1,s[a[p]]<u?l=p+1:o=p;u<s[a[l]]&&(l>0&&(n[e]=a[l-1]),a[l]=e)}}for(l=a.length,o=a[l-1];l-- >0;)a[l]=o,o=n[o];return a}function Hc(s){const n=s.subTree.component;if(n)return n.asyncDep&&!n.asyncResolved?n:Hc(n)}function ro(s){if(s)for(let n=0;n<s.length;n++)s[n].flags|=8}const Nu=Symbol.for("v-scx"),Bu=()=>on(Nu);function zu(s,n){return Tl(s,null,n)}function Bs(s,n,a){return Tl(s,n,a)}function Tl(s,n,a=ks){const{immediate:e,deep:t,flush:l,once:o}=a,p=Gs({},a),c=n&&e||!n&&l!=="post";let u;if(Ta){if(l==="sync"){const h=Bu();u=h.__watcherHandles||(h.__watcherHandles=[])}else if(!c){const h=()=>{};return h.stop=An,h.resume=An,h.pause=An,h}}const r=Ys;p.call=(h,j,m)=>wn(h,r,j,m);let i=!1;l==="post"?p.scheduler=h=>{an(h,r&&r.suspense)}:l!=="sync"&&(i=!0,p.scheduler=(h,j)=>{j?h():yl(h)}),p.augmentJob=h=>{n&&(h.flags|=4),i&&(h.flags|=2,r&&(h.id=r.uid,h.i=r))};const d=nu(s,n,p);return Ta&&(u?u.push(d):c&&d()),d}function qu(s,n,a){const e=this.proxy,t=Ms(s)?s.includes(".")?Vc(e,s):()=>e[s]:s.bind(e,e);let l;ds(n)?l=n:(l=n.handler,a=n);const o=je(this),p=Tl(t,l.bind(e),a);return o(),p}function Vc(s,n){const a=n.split(".");return()=>{let e=s;for(let t=0;t<a.length&&e;t++)e=e[a[t]];return e}}const Hu=(s,n)=>n==="modelValue"||n==="model-value"?s.modelModifiers:s[`${n}Modifiers`]||s[`${yn(n)}Modifiers`]||s[`${la(n)}Modifiers`];function Vu(s,n,...a){if(s.isUnmounted)return;const e=s.vnode.props||ks;let t=a;const l=n.startsWith("update:"),o=l&&Hu(e,n.slice(7));o&&(o.trim&&(t=a.map(r=>Ms(r)?r.trim():r)),o.number&&(t=a.map(rl)));let p,c=e[p=dt(n)]||e[p=dt(yn(n))];!c&&l&&(c=e[p=dt(la(n))]),c&&wn(c,s,6,t);const u=e[p+"Once"];if(u){if(!s.emitted)s.emitted={};else if(s.emitted[p])return;s.emitted[p]=!0,wn(u,s,6,t)}}const Gu=new WeakMap;function Gc(s,n,a=!1){const e=a?Gu:n.emitsCache,t=e.get(s);if(t!==void 0)return t;const l=s.emits;let o={},p=!1;if(!ds(s)){const c=u=>{const r=Gc(u,n,!0);r&&(p=!0,Gs(o,r))};!a&&n.mixins.length&&n.mixins.forEach(c),s.extends&&c(s.extends),s.mixins&&s.mixins.forEach(c)}return!l&&!p?(Rs(s)&&e.set(s,null),null):(is(l)?l.forEach(c=>o[c]=null):Gs(o,l),Rs(s)&&e.set(s,o),o)}function et(s,n){return!s||!Xe(n)?!1:(n=n.slice(2).replace(/Once$/,""),xs(s,n[0].toLowerCase()+n.slice(1))||xs(s,la(n))||xs(s,n))}function io(s){const{type:n,vnode:a,proxy:e,withProxy:t,propsOptions:[l],slots:o,attrs:p,emit:c,render:u,renderCache:r,props:i,data:d,setupState:h,ctx:j,inheritAttrs:m}=s,P=Be(s);let C,b;try{if(a.shapeFlag&4){const f=t||e,k=f;C=Pn(u.call(k,f,r,i,h,d,j)),b=p}else{const f=n;C=Pn(f.length>1?f(i,{attrs:p,slots:o,emit:c}):f(i,null)),b=n.props?p:Uu(p)}}catch(f){Ya.length=0,ue(f,s,1),C=bs(ln)}let x=C;if(b&&m!==!1){const f=Object.keys(b),{shapeFlag:k}=x;f.length&&k&7&&(l&&f.some(ol)&&(b=Wu(b,l)),x=aa(x,b,!1,!0))}return a.dirs&&(x=aa(x,null,!1,!0),x.dirs=x.dirs?x.dirs.concat(a.dirs):a.dirs),a.transition&&ae(x,a.transition),C=x,Be(P),C}const Uu=s=>{let n;for(const a in s)(a==="class"||a==="style"||Xe(a))&&((n||(n={}))[a]=s[a]);return n},Wu=(s,n)=>{const a={};for(const e in s)(!ol(e)||!(e.slice(9)in n))&&(a[e]=s[e]);return a};function Ku(s,n,a){const{props:e,children:t,component:l}=s,{props:o,children:p,patchFlag:c}=n,u=l.emitsOptions;if(n.dirs||n.transition)return!0;if(a&&c>=0){if(c&1024)return!0;if(c&16)return e?uo(e,o,u):!!o;if(c&8){const r=n.dynamicProps;for(let i=0;i<r.length;i++){const d=r[i];if(o[d]!==e[d]&&!et(u,d))return!0}}}else return(t||p)&&(!p||!p.$stable)?!0:e===o?!1:e?o?uo(e,o,u):!0:!!o;return!1}function uo(s,n,a){const e=Object.keys(n);if(e.length!==Object.keys(s).length)return!0;for(let t=0;t<e.length;t++){const l=e[t];if(n[l]!==s[l]&&!et(a,l))return!0}return!1}function Yu({vnode:s,parent:n},a){for(;n;){const e=n.subTree;if(e.suspense&&e.suspense.activeBranch===s&&(e.el=s.el),e===s)(s=n.vnode).el=a,n=n.parent;else break}}const Uc=s=>s.__isSuspense;function Xu(s,n){n&&n.pendingBranch?is(s)?n.effects.push(...s):n.effects.push(s):tu(s)}const js=Symbol.for("v-fgt"),me=Symbol.for("v-txt"),ln=Symbol.for("v-cmt"),Ee=Symbol.for("v-stc"),Ya=[];let rn=null;function E(s=!1){Ya.push(rn=s?null:[])}function Ju(){Ya.pop(),rn=Ya[Ya.length-1]||null}let ee=1;function He(s,n=!1){ee+=s,s<0&&rn&&n&&(rn.hasOnce=!0)}function Wc(s){return s.dynamicChildren=ee>0?rn||Ca:null,Ju(),ee>0&&rn&&rn.push(s),s}function L(s,n,a,e,t,l){return Wc(v(s,n,a,e,t,l,!0))}function ma(s,n,a,e,t){return Wc(bs(s,n,a,e,t,!0))}function Ve(s){return s?s.__v_isVNode===!0:!1}function da(s,n){return s.type===n.type&&s.key===n.key}const Kc=({key:s})=>s??null,Me=({ref:s,ref_key:n,ref_for:a})=>(typeof s=="number"&&(s=""+s),s!=null?Ms(s)||Es(s)||ds(s)?{i:cn,r:s,k:n,f:!!a}:s:null);function v(s,n=null,a=null,e=0,t=null,l=s===js?0:1,o=!1,p=!1){const c={__v_isVNode:!0,__v_skip:!0,type:s,props:n,key:n&&Kc(n),ref:n&&Me(n),scopeId:fc,slotScopeIds:null,children:a,component:null,suspense:null,ssContent:null,ssFallback:null,dirs:null,transition:null,el:null,anchor:null,target:null,targetStart:null,targetAnchor:null,staticCount:0,shapeFlag:l,patchFlag:e,dynamicProps:t,dynamicChildren:null,appContext:null,ctx:cn};return p?(Rl(c,a),l&128&&s.normalize(c)):a&&(c.shapeFlag|=Ms(a)?8:16),ee>0&&!o&&rn&&(c.patchFlag>0||l&6)&&c.patchFlag!==32&&rn.push(c),c}const bs=Qu;function Qu(s,n=null,a=null,e=0,t=null,l=!1){if((!s||s===Rc)&&(s=ln),Ve(s)){const p=aa(s,n,!0);return a&&Rl(p,a),ee>0&&!l&&rn&&(p.shapeFlag&6?rn[rn.indexOf(s)]=p:rn.push(p)),p.patchFlag=-2,p}if(rh(s)&&(s=s.__vccOpts),n){n=Zu(n);let{class:p,style:c}=n;p&&!Ms(p)&&(n.class=Ns(p)),Rs(c)&&(fl(c)&&!is(c)&&(c=Gs({},c)),n.style=oa(c))}const o=Ms(s)?1:Uc(s)?128:_c(s)?64:Rs(s)?4:ds(s)?2:0;return v(s,n,a,e,t,o,l,!0)}function Zu(s){return s?fl(s)||$c(s)?Gs({},s):s:null}function aa(s,n,a=!1,e=!1){const{props:t,ref:l,patchFlag:o,children:p,transition:c}=s,u=n?sh(t||{},n):t,r={__v_isVNode:!0,__v_skip:!0,type:s.type,props:u,key:u&&Kc(u),ref:n&&n.ref?a&&l?is(l)?l.concat(Me(n)):[l,Me(n)]:Me(n):l,scopeId:s.scopeId,slotScopeIds:s.slotScopeIds,children:p,target:s.target,targetStart:s.targetStart,targetAnchor:s.targetAnchor,staticCount:s.staticCount,shapeFlag:s.shapeFlag,patchFlag:n&&s.type!==js?o===-1?16:o|16:o,dynamicProps:s.dynamicProps,dynamicChildren:s.dynamicChildren,appContext:s.appContext,dirs:s.dirs,transition:c,component:s.component,suspense:s.suspense,ssContent:s.ssContent&&aa(s.ssContent),ssFallback:s.ssFallback&&aa(s.ssFallback),placeholder:s.placeholder,el:s.el,anchor:s.anchor,ctx:s.ctx,ce:s.ce};return c&&e&&ae(r,c.clone(r)),r}function Ss(s=" ",n=0){return bs(me,null,s,n)}function Ce(s,n){const a=bs(Ee,null,s);return a.staticCount=n,a}function rs(s="",n=!1){return n?(E(),ma(ln,null,s)):bs(ln,null,s)}function Pn(s){return s==null||typeof s=="boolean"?bs(ln):is(s)?bs(js,null,s.slice()):Ve(s)?Xn(s):bs(me,null,String(s))}function Xn(s){return s.el===null&&s.patchFlag!==-1||s.memo?s:aa(s)}function Rl(s,n){let a=0;const{shapeFlag:e}=s;if(n==null)n=null;else if(is(n))a=16;else if(typeof n=="object")if(e&65){const t=n.default;t&&(t._c&&(t._d=!1),Rl(s,t()),t._c&&(t._d=!0));return}else{a=32;const t=n._;!t&&!$c(n)?n._ctx=cn:t===3&&cn&&(cn.slots._===1?n._=1:(n._=2,s.patchFlag|=1024))}else ds(n)?(n={default:n,_ctx:cn},a=32):(n=String(n),e&64?(a=16,n=[Ss(n)]):a=8);s.children=n,s.shapeFlag|=a}function sh(...s){const n={};for(let a=0;a<s.length;a++){const e=s[a];for(const t in e)if(t==="class")n.class!==e.class&&(n.class=Ns([n.class,e.class]));else if(t==="style")n.style=oa([n.style,e.style]);else if(Xe(t)){const l=n[t],o=e[t];o&&l!==o&&!(is(l)&&l.includes(o))&&(n[t]=l?[].concat(l,o):o)}else t!==""&&(n[t]=e[t])}return n}function kn(s,n,a,e=null){wn(s,n,7,[a,e])}const nh=Dc();let ah=0;function eh(s,n,a){const e=s.type,t=(n?n.appContext:s.appContext)||nh,l={uid:ah++,vnode:s,type:e,parent:n,appContext:t,root:null,next:null,subTree:null,effect:null,update:null,job:null,scope:new Kp(!0),render:null,proxy:null,exposed:null,exposeProxy:null,withProxy:null,provides:n?n.provides:Object.create(t.provides),ids:n?n.ids:["",0,0],accessCache:null,renderCache:[],components:null,directives:null,propsOptions:Nc(e,t),emitsOptions:Gc(e,t),emit:null,emitted:null,propsDefaults:ks,inheritAttrs:e.inheritAttrs,ctx:ks,data:ks,props:ks,attrs:ks,slots:ks,refs:ks,setupState:ks,setupContext:null,suspense:a,suspenseId:a?a.pendingId:0,asyncDep:null,asyncResolved:!1,isMounted:!1,isUnmounted:!1,isDeactivated:!1,bc:null,c:null,bm:null,m:null,bu:null,u:null,um:null,bum:null,da:null,a:null,rtg:null,rtc:null,ec:null,sp:null};return l.ctx={_:l},l.root=n?n.root:l,l.emit=Vu.bind(null,l),s.ce&&s.ce(l),l}let Ys=null;const Gn=()=>Ys||cn;let Ge,Gt;{const s=st(),n=(a,e)=>{let t;return(t=s[a])||(t=s[a]=[]),t.push(e),l=>{t.length>1?t.forEach(o=>o(l)):t[0](l)}};Ge=n("__VUE_INSTANCE_SETTERS__",a=>Ys=a),Gt=n("__VUE_SSR_SETTERS__",a=>Ta=a)}const je=s=>{const n=Ys;return Ge(s),s.scope.on(),()=>{s.scope.off(),Ge(n)}},ho=()=>{Ys&&Ys.scope.off(),Ge(null)};function Yc(s){return s.vnode.shapeFlag&4}let Ta=!1;function th(s,n=!1,a=!1){n&&Gt(n);const{props:e,children:t}=s.vnode,l=Yc(s);Tu(s,e,l,n),Du(s,t,a||n);const o=l?lh(s,n):void 0;return n&&Gt(!1),o}function lh(s,n){const a=s.type;s.accessCache=Object.create(null),s.proxy=new Proxy(s.ctx,vu);const{setup:e}=a;if(e){zn();const t=s.setupContext=e.length>1?ph(s):null,l=je(s),o=ie(e,s,0,[s.props,t]),p=zp(o);if(qn(),l(),(p||s.sp)&&!Wa(s)&&vl(s),p){if(o.then(ho,ho),n)return o.then(c=>{mo(s,c)}).catch(c=>{ue(c,s,0)});s.asyncDep=o}else mo(s,o)}else Xc(s)}function mo(s,n,a){ds(n)?s.type.__ssrInlineRender?s.ssrRender=n:s.render=n:Rs(n)&&(s.setupState=hc(n)),Xc(s)}function Xc(s,n,a){const e=s.type;s.render||(s.render=e.render||An);{const t=je(s);zn();try{wu(s)}finally{qn(),t()}}}const oh={get(s,n){return Qs(s,"get",""),s[n]}};function ph(s){const n=a=>{s.exposed=a||{}};return{attrs:new Proxy(s.attrs,oh),slots:s.slots,emit:s.emit,expose:n}}function tt(s){return s.exposed?s.exposeProxy||(s.exposeProxy=new Proxy(hc(bl(s.exposed)),{get(n,a){if(a in n)return n[a];if(a in Ka)return Ka[a](s)},has(n,a){return a in n||a in Ka}})):s.proxy}function ch(s,n=!0){return ds(s)?s.displayName||s.name:s.name||n&&s.__name}function rh(s){return ds(s)&&"__vccOpts"in s}const ss=(s,n)=>Zi(s,n,Ta);function ge(s,n,a){try{He(-1);const e=arguments.length;return e===2?Rs(n)&&!is(n)?Ve(n)?bs(s,null,[n]):bs(s,n):bs(s,null,n):(e>3?a=Array.prototype.slice.call(arguments,2):e===3&&Ve(a)&&(a=[a]),bs(s,n,a))}finally{He(1)}}const ih="3.5.24";/**
* @vue/runtime-dom v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/let Ut;const jo=typeof window<"u"&&window.trustedTypes;if(jo)try{Ut=jo.createPolicy("vue",{createHTML:s=>s})}catch{}const Jc=Ut?s=>Ut.createHTML(s):s=>s,uh="http://www.w3.org/2000/svg",hh="http://www.w3.org/1998/Math/MathML",Dn=typeof document<"u"?document:null,go=Dn&&Dn.createElement("template"),dh={insert:(s,n,a)=>{n.insertBefore(s,a||null)},remove:s=>{const n=s.parentNode;n&&n.removeChild(s)},createElement:(s,n,a,e)=>{const t=n==="svg"?Dn.createElementNS(uh,s):n==="mathml"?Dn.createElementNS(hh,s):a?Dn.createElement(s,{is:a}):Dn.createElement(s);return s==="select"&&e&&e.multiple!=null&&t.setAttribute("multiple",e.multiple),t},createText:s=>Dn.createTextNode(s),createComment:s=>Dn.createComment(s),setText:(s,n)=>{s.nodeValue=n},setElementText:(s,n)=>{s.textContent=n},parentNode:s=>s.parentNode,nextSibling:s=>s.nextSibling,querySelector:s=>Dn.querySelector(s),setScopeId(s,n){s.setAttribute(n,"")},insertStaticContent(s,n,a,e,t,l){const o=a?a.previousSibling:n.lastChild;if(t&&(t===l||t.nextSibling))for(;n.insertBefore(t.cloneNode(!0),a),!(t===l||!(t=t.nextSibling)););else{go.innerHTML=Jc(e==="svg"?`<svg>${s}</svg>`:e==="mathml"?`<math>${s}</math>`:s);const p=go.content;if(e==="svg"||e==="mathml"){const c=p.firstChild;for(;c.firstChild;)p.appendChild(c.firstChild);p.removeChild(c)}n.insertBefore(p,a)}return[o?o.nextSibling:n.firstChild,a?a.previousSibling:n.lastChild]}},Un="transition",Fa="animation",te=Symbol("_vtc"),Qc={name:String,type:String,css:{type:Boolean,default:!0},duration:[String,Number,Object],enterFromClass:String,enterActiveClass:String,enterToClass:String,appearFromClass:String,appearActiveClass:String,appearToClass:String,leaveFromClass:String,leaveActiveClass:String,leaveToClass:String},mh=Gs({},wc,Qc),jh=s=>(s.displayName="Transition",s.props=mh,s),gh=jh((s,{slots:n})=>ge(cu,fh(s),n)),ia=(s,n=[])=>{is(s)?s.forEach(a=>a(...n)):s&&s(...n)},fo=s=>s?is(s)?s.some(n=>n.length>1):s.length>1:!1;function fh(s){const n={};for(const B in s)B in Qc||(n[B]=s[B]);if(s.css===!1)return n;const{name:a="v",type:e,duration:t,enterFromClass:l=`${a}-enter-from`,enterActiveClass:o=`${a}-enter-active`,enterToClass:p=`${a}-enter-to`,appearFromClass:c=l,appearActiveClass:u=o,appearToClass:r=p,leaveFromClass:i=`${a}-leave-from`,leaveActiveClass:d=`${a}-leave-active`,leaveToClass:h=`${a}-leave-to`}=s,j=bh(t),m=j&&j[0],P=j&&j[1],{onBeforeEnter:C,onEnter:b,onEnterCancelled:x,onLeave:f,onLeaveCancelled:k,onBeforeAppear:w=C,onAppear:y=b,onAppearCancelled:M=x}=n,S=(B,O,V,us)=>{B._enterCancelled=us,ua(B,O?r:p),ua(B,O?u:o),V&&V()},F=(B,O)=>{B._isLeaving=!1,ua(B,i),ua(B,h),ua(B,d),O&&O()},Z=B=>(O,V)=>{const us=B?y:b,os=()=>S(O,B,V);ia(us,[O,os]),bo(()=>{ua(O,B?c:l),En(O,B?r:p),fo(us)||_o(O,e,m,os)})};return Gs(n,{onBeforeEnter(B){ia(C,[B]),En(B,l),En(B,o)},onBeforeAppear(B){ia(w,[B]),En(B,c),En(B,u)},onEnter:Z(!1),onAppear:Z(!0),onLeave(B,O){B._isLeaving=!0;const V=()=>F(B,O);En(B,i),B._enterCancelled?(En(B,d),wo(B)):(wo(B),En(B,d)),bo(()=>{B._isLeaving&&(ua(B,i),En(B,h),fo(f)||_o(B,e,P,V))}),ia(f,[B,V])},onEnterCancelled(B){S(B,!1,void 0,!0),ia(x,[B])},onAppearCancelled(B){S(B,!0,void 0,!0),ia(M,[B])},onLeaveCancelled(B){F(B),ia(k,[B])}})}function bh(s){if(s==null)return null;if(Rs(s))return[vt(s.enter),vt(s.leave)];{const n=vt(s);return[n,n]}}function vt(s){return _i(s)}function En(s,n){n.split(/\s+/).forEach(a=>a&&s.classList.add(a)),(s[te]||(s[te]=new Set)).add(n)}function ua(s,n){n.split(/\s+/).forEach(e=>e&&s.classList.remove(e));const a=s[te];a&&(a.delete(n),a.size||(s[te]=void 0))}function bo(s){requestAnimationFrame(()=>{requestAnimationFrame(s)})}let _h=0;function _o(s,n,a,e){const t=s._endId=++_h,l=()=>{t===s._endId&&e()};if(a!=null)return setTimeout(l,a);const{type:o,timeout:p,propCount:c}=yh(s,n);if(!o)return e();const u=o+"end";let r=0;const i=()=>{s.removeEventListener(u,d),l()},d=h=>{h.target===s&&++r>=c&&i()};setTimeout(()=>{r<c&&i()},p+1),s.addEventListener(u,d)}function yh(s,n){const a=window.getComputedStyle(s),e=j=>(a[j]||"").split(", "),t=e(`${Un}Delay`),l=e(`${Un}Duration`),o=yo(t,l),p=e(`${Fa}Delay`),c=e(`${Fa}Duration`),u=yo(p,c);let r=null,i=0,d=0;n===Un?o>0&&(r=Un,i=o,d=l.length):n===Fa?u>0&&(r=Fa,i=u,d=c.length):(i=Math.max(o,u),r=i>0?o>u?Un:Fa:null,d=r?r===Un?l.length:c.length:0);const h=r===Un&&/\b(?:transform|all)(?:,|$)/.test(e(`${Un}Property`).toString());return{type:r,timeout:i,propCount:d,hasTransform:h}}function yo(s,n){for(;s.length<n.length;)s=s.concat(s);return Math.max(...n.map((a,e)=>vo(a)+vo(s[e])))}function vo(s){return s==="auto"?0:Number(s.slice(0,-1).replace(",","."))*1e3}function wo(s){return(s?s.ownerDocument:document).body.offsetHeight}function vh(s,n,a){const e=s[te];e&&(n=(n?[n,...e]:[...e]).join(" ")),n==null?s.removeAttribute("class"):a?s.setAttribute("class",n):s.className=n}const Ue=Symbol("_vod"),Zc=Symbol("_vsh"),sr={name:"show",beforeMount(s,{value:n},{transition:a}){s[Ue]=s.style.display==="none"?"":s.style.display,a&&n?a.beforeEnter(s):$a(s,n)},mounted(s,{value:n},{transition:a}){a&&n&&a.enter(s)},updated(s,{value:n,oldValue:a},{transition:e}){!n!=!a&&(e?n?(e.beforeEnter(s),$a(s,!0),e.enter(s)):e.leave(s,()=>{$a(s,!1)}):$a(s,n))},beforeUnmount(s,{value:n}){$a(s,n)}};function $a(s,n){s.style.display=n?s[Ue]:"none",s[Zc]=!n}const wh=Symbol(""),xh=/(?:^|;)\s*display\s*:/;function Ch(s,n,a){const e=s.style,t=Ms(a);let l=!1;if(a&&!t){if(n)if(Ms(n))for(const o of n.split(";")){const p=o.slice(0,o.indexOf(":")).trim();a[p]==null&&De(e,p,"")}else for(const o in n)a[o]==null&&De(e,o,"");for(const o in a)o==="display"&&(l=!0),De(e,o,a[o])}else if(t){if(n!==a){const o=e[wh];o&&(a+=";"+o),e.cssText=a,l=xh.test(a)}}else n&&s.removeAttribute("style");Ue in s&&(s[Ue]=l?e.display:"",s[Zc]&&(e.display="none"))}const xo=/\s*!important$/;function De(s,n,a){if(is(a))a.forEach(e=>De(s,n,e));else if(a==null&&(a=""),n.startsWith("--"))s.setProperty(n,a);else{const e=kh(s,n);xo.test(a)?s.setProperty(la(e),a.replace(xo,""),"important"):s[e]=a}}const Co=["Webkit","Moz","ms"],wt={};function kh(s,n){const a=wt[n];if(a)return a;let e=yn(n);if(e!=="filter"&&e in s)return wt[n]=e;e=Ze(e);for(let t=0;t<Co.length;t++){const l=Co[t]+e;if(l in s)return wt[n]=l}return n}const ko="http://www.w3.org/1999/xlink";function So(s,n,a,e,t,l=ki(n)){e&&n.startsWith("xlink:")?a==null?s.removeAttributeNS(ko,n.slice(6,n.length)):s.setAttributeNS(ko,n,a):a==null||l&&!Gp(a)?s.removeAttribute(n):s.setAttribute(n,l?"":ta(a)?String(a):a)}function Po(s,n,a,e,t){if(n==="innerHTML"||n==="textContent"){a!=null&&(s[n]=n==="innerHTML"?Jc(a):a);return}const l=s.tagName;if(n==="value"&&l!=="PROGRESS"&&!l.includes("-")){const p=l==="OPTION"?s.getAttribute("value")||"":s.value,c=a==null?s.type==="checkbox"?"on":"":String(a);(p!==c||!("_value"in s))&&(s.value=c),a==null&&s.removeAttribute(n),s._value=a;return}let o=!1;if(a===""||a==null){const p=typeof s[n];p==="boolean"?a=Gp(a):a==null&&p==="string"?(a="",o=!0):p==="number"&&(a=0,o=!0)}try{s[n]=a}catch{}o&&s.removeAttribute(t||n)}function wa(s,n,a,e){s.addEventListener(n,a,e)}function Sh(s,n,a,e){s.removeEventListener(n,a,e)}const Ao=Symbol("_vei");function Ph(s,n,a,e,t=null){const l=s[Ao]||(s[Ao]={}),o=l[n];if(e&&o)o.value=e;else{const[p,c]=Ah(n);if(e){const u=l[n]=Eh(e,t);wa(s,p,u,c)}else o&&(Sh(s,p,o,c),l[n]=void 0)}}const To=/(?:Once|Passive|Capture)$/;function Ah(s){let n;if(To.test(s)){n={};let e;for(;e=s.match(To);)s=s.slice(0,s.length-e[0].length),n[e[0].toLowerCase()]=!0}return[s[2]===":"?s.slice(3):la(s.slice(2)),n]}let xt=0;const Th=Promise.resolve(),Rh=()=>xt||(Th.then(()=>xt=0),xt=Date.now());function Eh(s,n){const a=e=>{if(!e._vts)e._vts=Date.now();else if(e._vts<=a.attached)return;wn(Mh(e,a.value),n,5,[e])};return a.value=s,a.attached=Rh(),a}function Mh(s,n){if(is(n)){const a=s.stopImmediatePropagation;return s.stopImmediatePropagation=()=>{a.call(s),s._stopped=!0},n.map(e=>t=>!t._stopped&&e&&e(t))}else return n}const Ro=s=>s.charCodeAt(0)===111&&s.charCodeAt(1)===110&&s.charCodeAt(2)>96&&s.charCodeAt(2)<123,Dh=(s,n,a,e,t,l)=>{const o=t==="svg";n==="class"?vh(s,e,o):n==="style"?Ch(s,a,e):Xe(n)?ol(n)||Ph(s,n,a,e,l):(n[0]==="."?(n=n.slice(1),!0):n[0]==="^"?(n=n.slice(1),!1):Lh(s,n,e,o))?(Po(s,n,e),!s.tagName.includes("-")&&(n==="value"||n==="checked"||n==="selected")&&So(s,n,e,o,l,n!=="value")):s._isVueCE&&(/[A-Z]/.test(n)||!Ms(e))?Po(s,yn(n),e,l,n):(n==="true-value"?s._trueValue=e:n==="false-value"&&(s._falseValue=e),So(s,n,e,o))};function Lh(s,n,a,e){if(e)return!!(n==="innerHTML"||n==="textContent"||n in s&&Ro(n)&&ds(a));if(n==="spellcheck"||n==="draggable"||n==="translate"||n==="autocorrect"||n==="sandbox"&&s.tagName==="IFRAME"||n==="form"||n==="list"&&s.tagName==="INPUT"||n==="type"&&s.tagName==="TEXTAREA")return!1;if(n==="width"||n==="height"){const t=s.tagName;if(t==="IMG"||t==="VIDEO"||t==="CANVAS"||t==="SOURCE")return!1}return Ro(n)&&Ms(a)?!1:n in s}const Eo=s=>{const n=s.props["onUpdate:modelValue"]||!1;return is(n)?a=>Pe(n,a):n};function Oh(s){s.target.composing=!0}function Mo(s){const n=s.target;n.composing&&(n.composing=!1,n.dispatchEvent(new Event("input")))}const Ct=Symbol("_assign");function Do(s,n,a){return n&&(s=s.trim()),a&&(s=rl(s)),s}const mv={created(s,{modifiers:{lazy:n,trim:a,number:e}},t){s[Ct]=Eo(t);const l=e||t.props&&t.props.type==="number";wa(s,n?"change":"input",o=>{o.target.composing||s[Ct](Do(s.value,a,l))}),(a||l)&&wa(s,"change",()=>{s.value=Do(s.value,a,l)}),n||(wa(s,"compositionstart",Oh),wa(s,"compositionend",Mo),wa(s,"change",Mo))},mounted(s,{value:n}){s.value=n??""},beforeUpdate(s,{value:n,oldValue:a,modifiers:{lazy:e,trim:t,number:l}},o){if(s[Ct]=Eo(o),s.composing)return;const p=(l||s.type==="number")&&!/^0\d/.test(s.value)?rl(s.value):s.value,c=n??"";p!==c&&(document.activeElement===s&&s.type!=="range"&&(e&&n===a||t&&s.value.trim()===c)||(s.value=c))}},Fh=["ctrl","shift","alt","meta"],$h={stop:s=>s.stopPropagation(),prevent:s=>s.preventDefault(),self:s=>s.target!==s.currentTarget,ctrl:s=>!s.ctrlKey,shift:s=>!s.shiftKey,alt:s=>!s.altKey,meta:s=>!s.metaKey,left:s=>"button"in s&&s.button!==0,middle:s=>"button"in s&&s.button!==1,right:s=>"button"in s&&s.button!==2,exact:(s,n)=>Fh.some(a=>s[`${a}Key`]&&!n.includes(a))},un=(s,n)=>{const a=s._withMods||(s._withMods={}),e=n.join(".");return a[e]||(a[e]=((t,...l)=>{for(let o=0;o<n.length;o++){const p=$h[n[o]];if(p&&p(t,n))return}return s(t,...l)}))},Ih={esc:"escape",space:" ",up:"arrow-up",left:"arrow-left",right:"arrow-right",down:"arrow-down",delete:"backspace"},Lo=(s,n)=>{const a=s._withKeys||(s._withKeys={}),e=n.join(".");return a[e]||(a[e]=(t=>{if(!("key"in t))return;const l=la(t.key);if(n.some(o=>o===l||Ih[o]===l))return s(t)}))},Nh=Gs({patchProp:Dh},dh);let Oo;function Bh(){return Oo||(Oo=Ou(Nh))}const zh=((...s)=>{const n=Bh().createApp(...s),{mount:a}=n;return n.mount=e=>{const t=Hh(e);if(!t)return;const l=n._component;!ds(l)&&!l.render&&!l.template&&(l.template=t.innerHTML),t.nodeType===1&&(t.textContent="");const o=a(t,!1,qh(t));return t instanceof Element&&(t.removeAttribute("v-cloak"),t.setAttribute("data-v-app","")),o},n});function qh(s){if(s instanceof SVGElement)return"svg";if(typeof MathMLElement=="function"&&s instanceof MathMLElement)return"mathml"}function Hh(s){return Ms(s)?document.querySelector(s):s}const Vh=(s,n)=>Es(n)?Wi(n):n,nr="usehead";function Gh(s){return{install(a){a.config.globalProperties.$unhead=s,a.config.globalProperties.$head=s,a.provide(nr,s)}}.install}function Uh(){if(Lc()){const s=on(nr);if(s)return s}throw new Error("useHead() was called without provide context, ensure you call it through the setup() function.")}function lt(s,n={}){const a=n.head||Uh();return a.ssr?a.push(s||{},n):Wh(a,s,n)}function Wh(s,n,a={}){const e=ns(!1);let t;return zu(()=>{const o=e.value?{}:Oe(n,Vh);t?t.patch(o):t=s.push(o,a)}),Gn()&&(Tn(()=>{t.dispose()}),Ac(()=>{e.value=!0}),Pc(()=>{e.value=!1})),t}const Kh={created(){let s=!1;const n=Gn();if(!n)return;const a=n.type;!a||!("head"in a)||(s=typeof a.head=="function"?()=>a.head.call(n.proxy):a.head,s&&lt(s))}};function Yh(s={}){const n=di({domOptions:{render:mi(()=>Ip(n),a=>setTimeout(a,0))},...s});return n.install=Gh(n),n}/*!
  * vue-router v4.5.1
  * (c) 2025 Eduardo San Martin Morote
  * @license MIT
  */const xa=typeof document<"u";function ar(s){return typeof s=="object"||"displayName"in s||"props"in s||"__vccOpts"in s}function Xh(s){return s.__esModule||s[Symbol.toStringTag]==="Module"||s.default&&ar(s.default)}const ws=Object.assign;function kt(s,n){const a={};for(const e in n){const t=n[e];a[e]=xn(t)?t.map(s):s(t)}return a}const Xa=()=>{},xn=Array.isArray,er=/#/g,Jh=/&/g,Qh=/\//g,Zh=/=/g,sd=/\?/g,tr=/\+/g,nd=/%5B/g,ad=/%5D/g,lr=/%5E/g,ed=/%60/g,or=/%7B/g,td=/%7C/g,pr=/%7D/g,ld=/%20/g;function El(s){return encodeURI(""+s).replace(td,"|").replace(nd,"[").replace(ad,"]")}function od(s){return El(s).replace(or,"{").replace(pr,"}").replace(lr,"^")}function Wt(s){return El(s).replace(tr,"%2B").replace(ld,"+").replace(er,"%23").replace(Jh,"%26").replace(ed,"`").replace(or,"{").replace(pr,"}").replace(lr,"^")}function pd(s){return Wt(s).replace(Zh,"%3D")}function cd(s){return El(s).replace(er,"%23").replace(sd,"%3F")}function rd(s){return s==null?"":cd(s).replace(Qh,"%2F")}function le(s){try{return decodeURIComponent(""+s)}catch{}return""+s}const id=/\/$/,ud=s=>s.replace(id,"");function St(s,n,a="/"){let e,t={},l="",o="";const p=n.indexOf("#");let c=n.indexOf("?");return p<c&&p>=0&&(c=-1),c>-1&&(e=n.slice(0,c),l=n.slice(c+1,p>-1?p:n.length),t=s(l)),p>-1&&(e=e||n.slice(0,p),o=n.slice(p,n.length)),e=jd(e??n,a),{fullPath:e+(l&&"?")+l+o,path:e,query:t,hash:le(o)}}function hd(s,n){const a=n.query?s(n.query):"";return n.path+(a&&"?")+a+(n.hash||"")}function Fo(s,n){return!n||!s.toLowerCase().startsWith(n.toLowerCase())?s:s.slice(n.length)||"/"}function dd(s,n,a){const e=n.matched.length-1,t=a.matched.length-1;return e>-1&&e===t&&Ra(n.matched[e],a.matched[t])&&cr(n.params,a.params)&&s(n.query)===s(a.query)&&n.hash===a.hash}function Ra(s,n){return(s.aliasOf||s)===(n.aliasOf||n)}function cr(s,n){if(Object.keys(s).length!==Object.keys(n).length)return!1;for(const a in s)if(!md(s[a],n[a]))return!1;return!0}function md(s,n){return xn(s)?$o(s,n):xn(n)?$o(n,s):s===n}function $o(s,n){return xn(n)?s.length===n.length&&s.every((a,e)=>a===n[e]):s.length===1&&s[0]===n}function jd(s,n){if(s.startsWith("/"))return s;if(!s)return n;const a=n.split("/"),e=s.split("/"),t=e[e.length-1];(t===".."||t===".")&&e.push("");let l=a.length-1,o,p;for(o=0;o<e.length;o++)if(p=e[o],p!==".")if(p==="..")l>1&&l--;else break;return a.slice(0,l).join("/")+"/"+e.slice(o).join("/")}const Wn={path:"/",name:void 0,params:{},query:{},hash:"",fullPath:"/",matched:[],meta:{},redirectedFrom:void 0};var oe;(function(s){s.pop="pop",s.push="push"})(oe||(oe={}));var Ja;(function(s){s.back="back",s.forward="forward",s.unknown=""})(Ja||(Ja={}));function gd(s){if(!s)if(xa){const n=document.querySelector("base");s=n&&n.getAttribute("href")||"/",s=s.replace(/^\w+:\/\/[^\/]+/,"")}else s="/";return s[0]!=="/"&&s[0]!=="#"&&(s="/"+s),ud(s)}const fd=/^[^#]+#/;function bd(s,n){return s.replace(fd,"#")+n}function _d(s,n){const a=document.documentElement.getBoundingClientRect(),e=s.getBoundingClientRect();return{behavior:n.behavior,left:e.left-a.left-(n.left||0),top:e.top-a.top-(n.top||0)}}const ot=()=>({left:window.scrollX,top:window.scrollY});function yd(s){let n;if("el"in s){const a=s.el,e=typeof a=="string"&&a.startsWith("#"),t=typeof a=="string"?e?document.getElementById(a.slice(1)):document.querySelector(a):a;if(!t)return;n=_d(t,s)}else n=s;"scrollBehavior"in document.documentElement.style?window.scrollTo(n):window.scrollTo(n.left!=null?n.left:window.scrollX,n.top!=null?n.top:window.scrollY)}function Io(s,n){return(history.state?history.state.position-n:-1)+s}const Kt=new Map;function vd(s,n){Kt.set(s,n)}function wd(s){const n=Kt.get(s);return Kt.delete(s),n}let xd=()=>location.protocol+"//"+location.host;function rr(s,n){const{pathname:a,search:e,hash:t}=n,l=s.indexOf("#");if(l>-1){let p=t.includes(s.slice(l))?s.slice(l).length:1,c=t.slice(p);return c[0]!=="/"&&(c="/"+c),Fo(c,"")}return Fo(a,s)+e+t}function Cd(s,n,a,e){let t=[],l=[],o=null;const p=({state:d})=>{const h=rr(s,location),j=a.value,m=n.value;let P=0;if(d){if(a.value=h,n.value=d,o&&o===j){o=null;return}P=m?d.position-m.position:0}else e(h);t.forEach(C=>{C(a.value,j,{delta:P,type:oe.pop,direction:P?P>0?Ja.forward:Ja.back:Ja.unknown})})};function c(){o=a.value}function u(d){t.push(d);const h=()=>{const j=t.indexOf(d);j>-1&&t.splice(j,1)};return l.push(h),h}function r(){const{history:d}=window;d.state&&d.replaceState(ws({},d.state,{scroll:ot()}),"")}function i(){for(const d of l)d();l=[],window.removeEventListener("popstate",p),window.removeEventListener("beforeunload",r)}return window.addEventListener("popstate",p),window.addEventListener("beforeunload",r,{passive:!0}),{pauseListeners:c,listen:u,destroy:i}}function No(s,n,a,e=!1,t=!1){return{back:s,current:n,forward:a,replaced:e,position:window.history.length,scroll:t?ot():null}}function kd(s){const{history:n,location:a}=window,e={value:rr(s,a)},t={value:n.state};t.value||l(e.value,{back:null,current:e.value,forward:null,position:n.length-1,replaced:!0,scroll:null},!0);function l(c,u,r){const i=s.indexOf("#"),d=i>-1?(a.host&&document.querySelector("base")?s:s.slice(i))+c:xd()+s+c;try{n[r?"replaceState":"pushState"](u,"",d),t.value=u}catch(h){console.error(h),a[r?"replace":"assign"](d)}}function o(c,u){const r=ws({},n.state,No(t.value.back,c,t.value.forward,!0),u,{position:t.value.position});l(c,r,!0),e.value=c}function p(c,u){const r=ws({},t.value,n.state,{forward:c,scroll:ot()});l(r.current,r,!0);const i=ws({},No(e.value,c,null),{position:r.position+1},u);l(c,i,!1),e.value=c}return{location:e,state:t,push:p,replace:o}}function Sd(s){s=gd(s);const n=kd(s),a=Cd(s,n.state,n.location,n.replace);function e(l,o=!0){o||a.pauseListeners(),history.go(l)}const t=ws({location:"",base:s,go:e,createHref:bd.bind(null,s)},n,a);return Object.defineProperty(t,"location",{enumerable:!0,get:()=>n.location.value}),Object.defineProperty(t,"state",{enumerable:!0,get:()=>n.state.value}),t}function Pd(s){return typeof s=="string"||s&&typeof s=="object"}function ir(s){return typeof s=="string"||typeof s=="symbol"}const ur=Symbol("");var Bo;(function(s){s[s.aborted=4]="aborted",s[s.cancelled=8]="cancelled",s[s.duplicated=16]="duplicated"})(Bo||(Bo={}));function Ea(s,n){return ws(new Error,{type:s,[ur]:!0},n)}function Mn(s,n){return s instanceof Error&&ur in s&&(n==null||!!(s.type&n))}const zo="[^/]+?",Ad={sensitive:!1,strict:!1,start:!0,end:!0},Td=/[.+*?^${}()[\]/\\]/g;function Rd(s,n){const a=ws({},Ad,n),e=[];let t=a.start?"^":"";const l=[];for(const u of s){const r=u.length?[]:[90];a.strict&&!u.length&&(t+="/");for(let i=0;i<u.length;i++){const d=u[i];let h=40+(a.sensitive?.25:0);if(d.type===0)i||(t+="/"),t+=d.value.replace(Td,"\\$&"),h+=40;else if(d.type===1){const{value:j,repeatable:m,optional:P,regexp:C}=d;l.push({name:j,repeatable:m,optional:P});const b=C||zo;if(b!==zo){h+=10;try{new RegExp(`(${b})`)}catch(f){throw new Error(`Invalid custom RegExp for param "${j}" (${b}): `+f.message)}}let x=m?`((?:${b})(?:/(?:${b}))*)`:`(${b})`;i||(x=P&&u.length<2?`(?:/${x})`:"/"+x),P&&(x+="?"),t+=x,h+=20,P&&(h+=-8),m&&(h+=-20),b===".*"&&(h+=-50)}r.push(h)}e.push(r)}if(a.strict&&a.end){const u=e.length-1;e[u][e[u].length-1]+=.7000000000000001}a.strict||(t+="/?"),a.end?t+="$":a.strict&&!t.endsWith("/")&&(t+="(?:/|$)");const o=new RegExp(t,a.sensitive?"":"i");function p(u){const r=u.match(o),i={};if(!r)return null;for(let d=1;d<r.length;d++){const h=r[d]||"",j=l[d-1];i[j.name]=h&&j.repeatable?h.split("/"):h}return i}function c(u){let r="",i=!1;for(const d of s){(!i||!r.endsWith("/"))&&(r+="/"),i=!1;for(const h of d)if(h.type===0)r+=h.value;else if(h.type===1){const{value:j,repeatable:m,optional:P}=h,C=j in u?u[j]:"";if(xn(C)&&!m)throw new Error(`Provided param "${j}" is an array but it is not repeatable (* or + modifiers)`);const b=xn(C)?C.join("/"):C;if(!b)if(P)d.length<2&&(r.endsWith("/")?r=r.slice(0,-1):i=!0);else throw new Error(`Missing required param "${j}"`);r+=b}}return r||"/"}return{re:o,score:e,keys:l,parse:p,stringify:c}}function Ed(s,n){let a=0;for(;a<s.length&&a<n.length;){const e=n[a]-s[a];if(e)return e;a++}return s.length<n.length?s.length===1&&s[0]===80?-1:1:s.length>n.length?n.length===1&&n[0]===80?1:-1:0}function hr(s,n){let a=0;const e=s.score,t=n.score;for(;a<e.length&&a<t.length;){const l=Ed(e[a],t[a]);if(l)return l;a++}if(Math.abs(t.length-e.length)===1){if(qo(e))return 1;if(qo(t))return-1}return t.length-e.length}function qo(s){const n=s[s.length-1];return s.length>0&&n[n.length-1]<0}const Md={type:0,value:""},Dd=/[a-zA-Z0-9_]/;function Ld(s){if(!s)return[[]];if(s==="/")return[[Md]];if(!s.startsWith("/"))throw new Error(`Invalid path "${s}"`);function n(h){throw new Error(`ERR (${a})/"${u}": ${h}`)}let a=0,e=a;const t=[];let l;function o(){l&&t.push(l),l=[]}let p=0,c,u="",r="";function i(){u&&(a===0?l.push({type:0,value:u}):a===1||a===2||a===3?(l.length>1&&(c==="*"||c==="+")&&n(`A repeatable param (${u}) must be alone in its segment. eg: '/:ids+.`),l.push({type:1,value:u,regexp:r,repeatable:c==="*"||c==="+",optional:c==="*"||c==="?"})):n("Invalid state to consume buffer"),u="")}function d(){u+=c}for(;p<s.length;){if(c=s[p++],c==="\\"&&a!==2){e=a,a=4;continue}switch(a){case 0:c==="/"?(u&&i(),o()):c===":"?(i(),a=1):d();break;case 4:d(),a=e;break;case 1:c==="("?a=2:Dd.test(c)?d():(i(),a=0,c!=="*"&&c!=="?"&&c!=="+"&&p--);break;case 2:c===")"?r[r.length-1]=="\\"?r=r.slice(0,-1)+c:a=3:r+=c;break;case 3:i(),a=0,c!=="*"&&c!=="?"&&c!=="+"&&p--,r="";break;default:n("Unknown state");break}}return a===2&&n(`Unfinished custom RegExp for param "${u}"`),i(),o(),t}function Od(s,n,a){const e=Rd(Ld(s.path),a),t=ws(e,{record:s,parent:n,children:[],alias:[]});return n&&!t.record.aliasOf==!n.record.aliasOf&&n.children.push(t),t}function Fd(s,n){const a=[],e=new Map;n=Uo({strict:!1,end:!0,sensitive:!1},n);function t(i){return e.get(i)}function l(i,d,h){const j=!h,m=Vo(i);m.aliasOf=h&&h.record;const P=Uo(n,i),C=[m];if("alias"in i){const f=typeof i.alias=="string"?[i.alias]:i.alias;for(const k of f)C.push(Vo(ws({},m,{components:h?h.record.components:m.components,path:k,aliasOf:h?h.record:m})))}let b,x;for(const f of C){const{path:k}=f;if(d&&k[0]!=="/"){const w=d.record.path,y=w[w.length-1]==="/"?"":"/";f.path=d.record.path+(k&&y+k)}if(b=Od(f,d,P),h?h.alias.push(b):(x=x||b,x!==b&&x.alias.push(b),j&&i.name&&!Go(b)&&o(i.name)),dr(b)&&c(b),m.children){const w=m.children;for(let y=0;y<w.length;y++)l(w[y],b,h&&h.children[y])}h=h||b}return x?()=>{o(x)}:Xa}function o(i){if(ir(i)){const d=e.get(i);d&&(e.delete(i),a.splice(a.indexOf(d),1),d.children.forEach(o),d.alias.forEach(o))}else{const d=a.indexOf(i);d>-1&&(a.splice(d,1),i.record.name&&e.delete(i.record.name),i.children.forEach(o),i.alias.forEach(o))}}function p(){return a}function c(i){const d=Nd(i,a);a.splice(d,0,i),i.record.name&&!Go(i)&&e.set(i.record.name,i)}function u(i,d){let h,j={},m,P;if("name"in i&&i.name){if(h=e.get(i.name),!h)throw Ea(1,{location:i});P=h.record.name,j=ws(Ho(d.params,h.keys.filter(x=>!x.optional).concat(h.parent?h.parent.keys.filter(x=>x.optional):[]).map(x=>x.name)),i.params&&Ho(i.params,h.keys.map(x=>x.name))),m=h.stringify(j)}else if(i.path!=null)m=i.path,h=a.find(x=>x.re.test(m)),h&&(j=h.parse(m),P=h.record.name);else{if(h=d.name?e.get(d.name):a.find(x=>x.re.test(d.path)),!h)throw Ea(1,{location:i,currentLocation:d});P=h.record.name,j=ws({},d.params,i.params),m=h.stringify(j)}const C=[];let b=h;for(;b;)C.unshift(b.record),b=b.parent;return{name:P,path:m,params:j,matched:C,meta:Id(C)}}s.forEach(i=>l(i));function r(){a.length=0,e.clear()}return{addRoute:l,resolve:u,removeRoute:o,clearRoutes:r,getRoutes:p,getRecordMatcher:t}}function Ho(s,n){const a={};for(const e of n)e in s&&(a[e]=s[e]);return a}function Vo(s){const n={path:s.path,redirect:s.redirect,name:s.name,meta:s.meta||{},aliasOf:s.aliasOf,beforeEnter:s.beforeEnter,props:$d(s),children:s.children||[],instances:{},leaveGuards:new Set,updateGuards:new Set,enterCallbacks:{},components:"components"in s?s.components||null:s.component&&{default:s.component}};return Object.defineProperty(n,"mods",{value:{}}),n}function $d(s){const n={},a=s.props||!1;if("component"in s)n.default=a;else for(const e in s.components)n[e]=typeof a=="object"?a[e]:a;return n}function Go(s){for(;s;){if(s.record.aliasOf)return!0;s=s.parent}return!1}function Id(s){return s.reduce((n,a)=>ws(n,a.meta),{})}function Uo(s,n){const a={};for(const e in s)a[e]=e in n?n[e]:s[e];return a}function Nd(s,n){let a=0,e=n.length;for(;a!==e;){const l=a+e>>1;hr(s,n[l])<0?e=l:a=l+1}const t=Bd(s);return t&&(e=n.lastIndexOf(t,e-1)),e}function Bd(s){let n=s;for(;n=n.parent;)if(dr(n)&&hr(s,n)===0)return n}function dr({record:s}){return!!(s.name||s.components&&Object.keys(s.components).length||s.redirect)}function zd(s){const n={};if(s===""||s==="?")return n;const e=(s[0]==="?"?s.slice(1):s).split("&");for(let t=0;t<e.length;++t){const l=e[t].replace(tr," "),o=l.indexOf("="),p=le(o<0?l:l.slice(0,o)),c=o<0?null:le(l.slice(o+1));if(p in n){let u=n[p];xn(u)||(u=n[p]=[u]),u.push(c)}else n[p]=c}return n}function Wo(s){let n="";for(let a in s){const e=s[a];if(a=pd(a),e==null){e!==void 0&&(n+=(n.length?"&":"")+a);continue}(xn(e)?e.map(l=>l&&Wt(l)):[e&&Wt(e)]).forEach(l=>{l!==void 0&&(n+=(n.length?"&":"")+a,l!=null&&(n+="="+l))})}return n}function qd(s){const n={};for(const a in s){const e=s[a];e!==void 0&&(n[a]=xn(e)?e.map(t=>t==null?null:""+t):e==null?e:""+e)}return n}const Hd=Symbol(""),Ko=Symbol(""),pt=Symbol(""),Ml=Symbol(""),Yt=Symbol("");function Ia(){let s=[];function n(e){return s.push(e),()=>{const t=s.indexOf(e);t>-1&&s.splice(t,1)}}function a(){s=[]}return{add:n,list:()=>s.slice(),reset:a}}function Jn(s,n,a,e,t,l=o=>o()){const o=e&&(e.enterCallbacks[t]=e.enterCallbacks[t]||[]);return()=>new Promise((p,c)=>{const u=d=>{d===!1?c(Ea(4,{from:a,to:n})):d instanceof Error?c(d):Pd(d)?c(Ea(2,{from:n,to:d})):(o&&e.enterCallbacks[t]===o&&typeof d=="function"&&o.push(d),p())},r=l(()=>s.call(e&&e.instances[t],n,a,u));let i=Promise.resolve(r);s.length<3&&(i=i.then(u)),i.catch(d=>c(d))})}function Pt(s,n,a,e,t=l=>l()){const l=[];for(const o of s)for(const p in o.components){let c=o.components[p];if(!(n!=="beforeRouteEnter"&&!o.instances[p]))if(ar(c)){const r=(c.__vccOpts||c)[n];r&&l.push(Jn(r,a,e,o,p,t))}else{let u=c();l.push(()=>u.then(r=>{if(!r)throw new Error(`Couldn't resolve component "${p}" at "${o.path}"`);const i=Xh(r)?r.default:r;o.mods[p]=r,o.components[p]=i;const h=(i.__vccOpts||i)[n];return h&&Jn(h,a,e,o,p,t)()}))}}return l}function Yo(s){const n=on(pt),a=on(Ml),e=ss(()=>{const c=ls(s.to);return n.resolve(c)}),t=ss(()=>{const{matched:c}=e.value,{length:u}=c,r=c[u-1],i=a.matched;if(!r||!i.length)return-1;const d=i.findIndex(Ra.bind(null,r));if(d>-1)return d;const h=Xo(c[u-2]);return u>1&&Xo(r)===h&&i[i.length-1].path!==h?i.findIndex(Ra.bind(null,c[u-2])):d}),l=ss(()=>t.value>-1&&Kd(a.params,e.value.params)),o=ss(()=>t.value>-1&&t.value===a.matched.length-1&&cr(a.params,e.value.params));function p(c={}){if(Wd(c)){const u=n[ls(s.replace)?"replace":"push"](ls(s.to)).catch(Xa);return s.viewTransition&&typeof document<"u"&&"startViewTransition"in document&&document.startViewTransition(()=>u),u}return Promise.resolve()}return{route:e,href:ss(()=>e.value.href),isActive:l,isExactActive:o,navigate:p}}function Vd(s){return s.length===1?s[0]:s}const Gd=zs({name:"RouterLink",compatConfig:{MODE:3},props:{to:{type:[String,Object],required:!0},replace:Boolean,activeClass:String,exactActiveClass:String,custom:Boolean,ariaCurrentValue:{type:String,default:"page"},viewTransition:Boolean},useLink:Yo,setup(s,{slots:n}){const a=re(Yo(s)),{options:e}=on(pt),t=ss(()=>({[Jo(s.activeClass,e.linkActiveClass,"router-link-active")]:a.isActive,[Jo(s.exactActiveClass,e.linkExactActiveClass,"router-link-exact-active")]:a.isExactActive}));return()=>{const l=n.default&&Vd(n.default(a));return s.custom?l:ge("a",{"aria-current":a.isExactActive?s.ariaCurrentValue:null,href:a.href,onClick:a.navigate,class:t.value},l)}}}),Ud=Gd;function Wd(s){if(!(s.metaKey||s.altKey||s.ctrlKey||s.shiftKey)&&!s.defaultPrevented&&!(s.button!==void 0&&s.button!==0)){if(s.currentTarget&&s.currentTarget.getAttribute){const n=s.currentTarget.getAttribute("target");if(/\b_blank\b/i.test(n))return}return s.preventDefault&&s.preventDefault(),!0}}function Kd(s,n){for(const a in n){const e=n[a],t=s[a];if(typeof e=="string"){if(e!==t)return!1}else if(!xn(t)||t.length!==e.length||e.some((l,o)=>l!==t[o]))return!1}return!0}function Xo(s){return s?s.aliasOf?s.aliasOf.path:s.path:""}const Jo=(s,n,a)=>s??n??a,Yd=zs({name:"RouterView",inheritAttrs:!1,props:{name:{type:String,default:"default"},route:Object},compatConfig:{MODE:3},setup(s,{attrs:n,slots:a}){const e=on(Yt),t=ss(()=>s.route||e.value),l=on(Ko,0),o=ss(()=>{let u=ls(l);const{matched:r}=t.value;let i;for(;(i=r[u])&&!i.components;)u++;return u}),p=ss(()=>t.value.matched[o.value]);Re(Ko,ss(()=>o.value+1)),Re(Hd,p),Re(Yt,t);const c=ns();return Bs(()=>[c.value,p.value,s.name],([u,r,i],[d,h,j])=>{r&&(r.instances[i]=u,h&&h!==r&&u&&u===d&&(r.leaveGuards.size||(r.leaveGuards=h.leaveGuards),r.updateGuards.size||(r.updateGuards=h.updateGuards))),u&&r&&(!h||!Ra(r,h)||!d)&&(r.enterCallbacks[i]||[]).forEach(m=>m(u))},{flush:"post"}),()=>{const u=t.value,r=s.name,i=p.value,d=i&&i.components[r];if(!d)return Qo(a.default,{Component:d,route:u});const h=i.props[r],j=h?h===!0?u.params:typeof h=="function"?h(u):h:null,P=ge(d,ws({},j,n,{onVnodeUnmounted:C=>{C.component.isUnmounted&&(i.instances[r]=null)},ref:c}));return Qo(a.default,{Component:P,route:u})||P}}});function Qo(s,n){if(!s)return null;const a=s(n);return a.length===1?a[0]:a}const Xd=Yd;function Jd(s){const n=Fd(s.routes,s),a=s.parseQuery||zd,e=s.stringifyQuery||Wo,t=s.history,l=Ia(),o=Ia(),p=Ia(),c=_l(Wn);let u=Wn;xa&&s.scrollBehavior&&"scrollRestoration"in history&&(history.scrollRestoration="manual");const r=kt.bind(null,R=>""+R),i=kt.bind(null,rd),d=kt.bind(null,le);function h(R,H){let q,X;return ir(R)?(q=n.getRecordMatcher(R),X=H):X=R,n.addRoute(X,q)}function j(R){const H=n.getRecordMatcher(R);H&&n.removeRoute(H)}function m(){return n.getRoutes().map(R=>R.record)}function P(R){return!!n.getRecordMatcher(R)}function C(R,H){if(H=ws({},H||c.value),typeof R=="string"){const A=St(a,R,H.path),$=n.resolve({path:A.path},H),N=t.createHref(A.fullPath);return ws(A,$,{params:d($.params),hash:le(A.hash),redirectedFrom:void 0,href:N})}let q;if(R.path!=null)q=ws({},R,{path:St(a,R.path,H.path).path});else{const A=ws({},R.params);for(const $ in A)A[$]==null&&delete A[$];q=ws({},R,{params:i(A)}),H.params=i(H.params)}const X=n.resolve(q,H),gs=R.hash||"";X.params=r(d(X.params));const g=hd(e,ws({},R,{hash:od(gs),path:X.path})),_=t.createHref(g);return ws({fullPath:g,hash:gs,query:e===Wo?qd(R.query):R.query||{}},X,{redirectedFrom:void 0,href:_})}function b(R){return typeof R=="string"?St(a,R,c.value.path):ws({},R)}function x(R,H){if(u!==R)return Ea(8,{from:H,to:R})}function f(R){return y(R)}function k(R){return f(ws(b(R),{replace:!0}))}function w(R){const H=R.matched[R.matched.length-1];if(H&&H.redirect){const{redirect:q}=H;let X=typeof q=="function"?q(R):q;return typeof X=="string"&&(X=X.includes("?")||X.includes("#")?X=b(X):{path:X},X.params={}),ws({query:R.query,hash:R.hash,params:X.path!=null?{}:R.params},X)}}function y(R,H){const q=u=C(R),X=c.value,gs=R.state,g=R.force,_=R.replace===!0,A=w(q);if(A)return y(ws(b(A),{state:typeof A=="object"?ws({},gs,A.state):gs,force:g,replace:_}),H||q);const $=q;$.redirectedFrom=H;let N;return!g&&dd(e,X,q)&&(N=Ea(16,{to:$,from:X}),Hs(X,X,!0,!1)),(N?Promise.resolve(N):F($,X)).catch(I=>Mn(I)?Mn(I,2)?I:$s(I):hs(I,$,X)).then(I=>{if(I){if(Mn(I,2))return y(ws({replace:_},b(I.to),{state:typeof I.to=="object"?ws({},gs,I.to.state):gs,force:g}),H||$)}else I=B($,X,!0,_,gs);return Z($,X,I),I})}function M(R,H){const q=x(R,H);return q?Promise.reject(q):Promise.resolve()}function S(R){const H=ms.values().next().value;return H&&typeof H.runWithContext=="function"?H.runWithContext(R):R()}function F(R,H){let q;const[X,gs,g]=Qd(R,H);q=Pt(X.reverse(),"beforeRouteLeave",R,H);for(const A of X)A.leaveGuards.forEach($=>{q.push(Jn($,R,H))});const _=M.bind(null,R,H);return q.push(_),fs(q).then(()=>{q=[];for(const A of l.list())q.push(Jn(A,R,H));return q.push(_),fs(q)}).then(()=>{q=Pt(gs,"beforeRouteUpdate",R,H);for(const A of gs)A.updateGuards.forEach($=>{q.push(Jn($,R,H))});return q.push(_),fs(q)}).then(()=>{q=[];for(const A of g)if(A.beforeEnter)if(xn(A.beforeEnter))for(const $ of A.beforeEnter)q.push(Jn($,R,H));else q.push(Jn(A.beforeEnter,R,H));return q.push(_),fs(q)}).then(()=>(R.matched.forEach(A=>A.enterCallbacks={}),q=Pt(g,"beforeRouteEnter",R,H,S),q.push(_),fs(q))).then(()=>{q=[];for(const A of o.list())q.push(Jn(A,R,H));return q.push(_),fs(q)}).catch(A=>Mn(A,8)?A:Promise.reject(A))}function Z(R,H,q){p.list().forEach(X=>S(()=>X(R,H,q)))}function B(R,H,q,X,gs){const g=x(R,H);if(g)return g;const _=H===Wn,A=xa?history.state:{};q&&(X||_?t.replace(R.fullPath,ws({scroll:_&&A&&A.scroll},gs)):t.push(R.fullPath,gs)),c.value=R,Hs(R,H,q,_),$s()}let O;function V(){O||(O=t.listen((R,H,q)=>{if(!_s.listening)return;const X=C(R),gs=w(X);if(gs){y(ws(gs,{replace:!0,force:!0}),X).catch(Xa);return}u=X;const g=c.value;xa&&vd(Io(g.fullPath,q.delta),ot()),F(X,g).catch(_=>Mn(_,12)?_:Mn(_,2)?(y(ws(b(_.to),{force:!0}),X).then(A=>{Mn(A,20)&&!q.delta&&q.type===oe.pop&&t.go(-1,!1)}).catch(Xa),Promise.reject()):(q.delta&&t.go(-q.delta,!1),hs(_,X,g))).then(_=>{_=_||B(X,g,!1),_&&(q.delta&&!Mn(_,8)?t.go(-q.delta,!1):q.type===oe.pop&&Mn(_,20)&&t.go(-1,!1)),Z(X,g,_)}).catch(Xa)}))}let us=Ia(),os=Ia(),Y;function hs(R,H,q){$s(R);const X=os.list();return X.length?X.forEach(gs=>gs(R,H,q)):console.error(R),Promise.reject(R)}function Fs(){return Y&&c.value!==Wn?Promise.resolve():new Promise((R,H)=>{us.add([R,H])})}function $s(R){return Y||(Y=!R,V(),us.list().forEach(([H,q])=>R?q(R):H()),us.reset()),R}function Hs(R,H,q,X){const{scrollBehavior:gs}=s;if(!xa||!gs)return Promise.resolve();const g=!q&&wd(Io(R.fullPath,0))||(X||!q)&&history.state&&history.state.scroll||null;return Ks().then(()=>gs(R,H,g)).then(_=>_&&yd(_)).catch(_=>hs(_,R,H))}const Vs=R=>t.go(R);let es;const ms=new Set,_s={currentRoute:c,listening:!0,addRoute:h,removeRoute:j,clearRoutes:n.clearRoutes,hasRoute:P,getRoutes:m,resolve:C,options:s,push:f,replace:k,go:Vs,back:()=>Vs(-1),forward:()=>Vs(1),beforeEach:l.add,beforeResolve:o.add,afterEach:p.add,onError:os.add,isReady:Fs,install(R){const H=this;R.component("RouterLink",Ud),R.component("RouterView",Xd),R.config.globalProperties.$router=H,Object.defineProperty(R.config.globalProperties,"$route",{enumerable:!0,get:()=>ls(c)}),xa&&!es&&c.value===Wn&&(es=!0,f(t.location).catch(gs=>{}));const q={};for(const gs in Wn)Object.defineProperty(q,gs,{get:()=>c.value[gs],enumerable:!0});R.provide(pt,H),R.provide(Ml,ic(q)),R.provide(Yt,c);const X=R.unmount;ms.add(R),R.unmount=function(){ms.delete(R),ms.size<1&&(u=Wn,O&&O(),O=null,c.value=Wn,es=!1,Y=!1),X()}}};function fs(R){return R.reduce((H,q)=>H.then(()=>S(q)),Promise.resolve())}return _s}function Qd(s,n){const a=[],e=[],t=[],l=Math.max(n.matched.length,s.matched.length);for(let o=0;o<l;o++){const p=n.matched[o];p&&(s.matched.find(u=>Ra(u,p))?e.push(p):a.push(p));const c=s.matched[o];c&&(n.matched.find(u=>Ra(u,c))||t.push(c))}return[a,e,t]}function fe(){return on(pt)}function Ma(s){return on(Ml)}function Zd(s){return document.readyState==="loading"?new Promise(n=>{document.addEventListener("DOMContentLoaded",()=>n(s))}):Promise.resolve(s)}const sm=zs({setup(s,{slots:n}){const a=ns(!1);return dn(()=>a.value=!0),()=>a.value?n.default&&n.default({}):n.placeholder&&n.placeholder({})}});function nm(s){try{return JSON.parse(s||"{}")}catch(n){return console.error("[SSG] On state deserialization -",n,s),{}}}function am(s,n,a,e){const{transformState:t,registerComponents:l=!0,useHead:o=!0,rootContainer:p="#app"}={};async function c(u){const r=zh(s);let i;o&&r.use(i=Yh());const d=Jd({history:Sd(n.base),...n}),{routes:h}=n;l&&r.component("ClientOnly",sm);const j=[],C={app:r,head:i,isClient:!0,router:d,routes:h,onSSRAppRendered:()=>{},triggerOnSSRAppRendered:()=>Promise.all(j.map(k=>k())),initialState:{},transformState:t,routePath:u};await Zd(),C.initialState=t?.(window.__INITIAL_STATE__||{})||nm(window.__INITIAL_STATE__),await a?.(C),r.use(d);let b,x=!0;d.beforeEach((k,w,y)=>{(x||b&&b===k.path)&&(x=!1,b=k.path,k.meta.state=C.initialState),y()});const f=C.initialState;return{...C,initialState:f}}return(async()=>{const{app:u,router:r}=await c();await r.isReady(),u.mount(p,!0)})(),c}/*!
 * pinia v3.0.3
 * (c) 2025 Eduardo San Martin Morote
 * @license MIT
 */let mr;const ct=s=>mr=s,jr=Symbol();function Xt(s){return s&&typeof s=="object"&&Object.prototype.toString.call(s)==="[object Object]"&&typeof s.toJSON!="function"}var Qa;(function(s){s.direct="direct",s.patchObject="patch object",s.patchFunction="patch function"})(Qa||(Qa={}));function em(){const s=il(!0),n=s.run(()=>ns({}));let a=[],e=[];const t=bl({install(l){ct(t),t._a=l,l.provide(jr,t),l.config.globalProperties.$pinia=t,e.forEach(o=>a.push(o)),e=[]},use(l){return this._a?a.push(l):e.push(l),this},_p:a,_a:null,_e:s,_s:new Map,state:n});return t}const gr=()=>{};function Zo(s,n,a,e=gr){s.push(n);const t=()=>{const l=s.indexOf(n);l>-1&&(s.splice(l,1),e())};return!a&&Yp()&&Si(t),t}function ya(s,...n){s.slice().forEach(a=>{a(...n)})}const tm=s=>s(),sp=Symbol(),At=Symbol();function Jt(s,n){s instanceof Map&&n instanceof Map?n.forEach((a,e)=>s.set(e,a)):s instanceof Set&&n instanceof Set&&n.forEach(s.add,s);for(const a in n){if(!n.hasOwnProperty(a))continue;const e=n[a],t=s[a];Xt(t)&&Xt(e)&&s.hasOwnProperty(a)&&!Es(e)&&!sa(e)?s[a]=Jt(t,e):s[a]=e}return s}const lm=Symbol();function om(s){return!Xt(s)||!Object.prototype.hasOwnProperty.call(s,lm)}const{assign:Kn}=Object;function pm(s){return!!(Es(s)&&s.effect)}function cm(s,n,a,e){const{state:t,actions:l,getters:o}=n,p=a.state.value[s];let c;function u(){p||(a.state.value[s]=t?t():{});const r=Yi(a.state.value[s]);return Kn(r,l,Object.keys(o||{}).reduce((i,d)=>(i[d]=bl(ss(()=>{ct(a);const h=a._s.get(s);return o[d].call(h,h)})),i),{}))}return c=fr(s,u,n,a,e,!0),c}function fr(s,n,a={},e,t,l){let o;const p=Kn({actions:{}},a),c={deep:!0};let u,r,i=[],d=[],h;const j=e.state.value[s];!l&&!j&&(e.state.value[s]={}),ns({});let m;function P(M){let S;u=r=!1,typeof M=="function"?(M(e.state.value[s]),S={type:Qa.patchFunction,storeId:s,events:h}):(Jt(e.state.value[s],M),S={type:Qa.patchObject,payload:M,storeId:s,events:h});const F=m=Symbol();Ks().then(()=>{m===F&&(u=!0)}),r=!0,ya(i,S,e.state.value[s])}const C=l?function(){const{state:S}=a,F=S?S():{};this.$patch(Z=>{Kn(Z,F)})}:gr;function b(){o.stop(),i=[],d=[],e._s.delete(s)}const x=(M,S="")=>{if(sp in M)return M[At]=S,M;const F=function(){ct(e);const Z=Array.from(arguments),B=[],O=[];function V(Y){B.push(Y)}function us(Y){O.push(Y)}ya(d,{args:Z,name:F[At],store:k,after:V,onError:us});let os;try{os=M.apply(this&&this.$id===s?this:k,Z)}catch(Y){throw ya(O,Y),Y}return os instanceof Promise?os.then(Y=>(ya(B,Y),Y)).catch(Y=>(ya(O,Y),Promise.reject(Y))):(ya(B,os),os)};return F[sp]=!0,F[At]=S,F},f={_p:e,$id:s,$onAction:Zo.bind(null,d),$patch:P,$reset:C,$subscribe(M,S={}){const F=Zo(i,M,S.detached,()=>Z()),Z=o.run(()=>Bs(()=>e.state.value[s],B=>{(S.flush==="sync"?r:u)&&M({storeId:s,type:Qa.direct,events:h},B)},Kn({},c,S)));return F},$dispose:b},k=re(f);e._s.set(s,k);const y=(e._a&&e._a.runWithContext||tm)(()=>e._e.run(()=>(o=il()).run(()=>n({action:x}))));for(const M in y){const S=y[M];if(Es(S)&&!pm(S)||sa(S))l||(j&&om(S)&&(Es(S)?S.value=j[M]:Jt(S,j[M])),e.state.value[s][M]=S);else if(typeof S=="function"){const F=x(S,M);y[M]=F,p.actions[M]=S}}return Kn(k,y),Kn(ys(k),y),Object.defineProperty(k,"$state",{get:()=>e.state.value[s],set:M=>{P(S=>{Kn(S,M)})}}),e._p.forEach(M=>{Kn(k,o.run(()=>M({store:k,app:e._a,pinia:e,options:p})))}),j&&l&&a.hydrate&&a.hydrate(k.$state,j),u=!0,r=!0,k}/*! #__NO_SIDE_EFFECTS__ */function rm(s,n,a){let e;const t=typeof n=="function";e=t?a:n;function l(o,p){const c=Lc();return o=o||(c?on(jr,null):null),o&&ct(o),o=mr,o._s.has(s)||(t?fr(s,n,e,o):cm(s,e,o)),o._s.get(s)}return l.$id=s,l}/*!
  * shared v12.0.0-alpha.3
  * (c) 2016-present kazuya kawaguchi and contributors
  * Released under the MIT License.
  */const np=typeof window<"u",pa=(s,n=!1)=>n?Symbol.for(s):Symbol(s),im=(s,n,a)=>um({l:s,k:n,s:a}),um=s=>JSON.stringify(s).replace(/\u2028/g,"\\u2028").replace(/\u2029/g,"\\u2029").replace(/\u0027/g,"\\u0027"),Zs=s=>typeof s=="number"&&isFinite(s),hm=s=>Dl(s)==="[object Date]",We=s=>Dl(s)==="[object RegExp]",rt=s=>vs(s)&&Object.keys(s).length===0,sn=Object.assign,dm=Object.create,Ps=(s=null)=>dm(s);let ap;const Pa=()=>ap||(ap=typeof globalThis<"u"?globalThis:typeof self<"u"?self:typeof window<"u"?window:typeof global<"u"?global:Ps());function ep(s){return s.replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;").replace(/'/g,"&apos;")}const mm=Object.prototype.hasOwnProperty;function fa(s,n){return mm.call(s,n)}const Is=Array.isArray,Ls=s=>typeof s=="function",as=s=>typeof s=="string",Os=s=>typeof s=="boolean",Cs=s=>s!==null&&typeof s=="object",jm=s=>Cs(s)&&Ls(s.then)&&Ls(s.catch),br=Object.prototype.toString,Dl=s=>br.call(s),vs=s=>Dl(s)==="[object Object]",gm=s=>s==null?"":Is(s)||vs(s)&&s.toString===br?JSON.stringify(s,null,2):String(s);function fm(s,n=""){return s.reduce((a,e,t)=>t===0?a+e:a+n+e,"")}function bm(s,n){typeof console<"u"&&(console.warn("[intlify] "+s),n&&console.warn(n.stack))}const ke=s=>!Cs(s)||Is(s);function Le(s,n){if(ke(s)||ke(n))throw new Error("Invalid value");const a=[{src:s,des:n}];for(;a.length;){const{src:e,des:t}=a.pop();Object.keys(e).forEach(l=>{l!=="__proto__"&&(Cs(e[l])&&!Cs(t[l])&&(t[l]=Array.isArray(e[l])?[]:Ps()),ke(t[l])||ke(e[l])?t[l]=e[l]:a.push({src:e[l],des:t[l]}))})}}/*!
  * message-compiler v12.0.0-alpha.3
  * (c) 2016-present kazuya kawaguchi and contributors
  * Released under the MIT License.
  */const _m=17;function _r(s,n,a={}){const{domain:e,messages:t,args:l}=a,o=s,p=new SyntaxError(String(o));return p.code=s,p.domain=e,p}/*!
  * core-base v12.0.0-alpha.3
  * (c) 2016-present kazuya kawaguchi and contributors
  * Released under the MIT License.
  */function ym(){typeof __INTLIFY_PROD_DEVTOOLS__!="boolean"&&(Pa().__INTLIFY_PROD_DEVTOOLS__=!1),typeof __INTLIFY_DROP_MESSAGE_COMPILER__!="boolean"&&(Pa().__INTLIFY_DROP_MESSAGE_COMPILER__=!1)}function Bn(s){return Cs(s)&&km(s)===0&&(fa(s,"b")||fa(s,"body"))}const vm=["b","body"],wm=["c","cases"],xm=["s","static"],Cm=["i","items"],yr=["t","type"];function km(s){return Tm(s,yr)}const Sm=["v","value"],Pm=["m","modifier"],Am=["k","key"];function Tm(s,n,a){for(let e=0;e<n.length;e++){const t=n[e];if(fa(s,t)&&s[t]!=null)return s[t]}return a}const Rm=[...vm,...wm,...xm,...Cm,...Am,...Pm,...Sm,...yr];let pe=null;function Em(s){pe=s}function Mm(s,n,a){pe&&pe.emit("i18n:init",{timestamp:Date.now(),i18n:s,version:n,meta:a})}const Dm=Lm("function:translate");function Lm(s){return n=>pe&&pe.emit(s,n)}const In={INVALID_ARGUMENT:_m,INVALID_DATE_ARGUMENT:18,INVALID_ISO_DATE_ARGUMENT:19,NOT_SUPPORT_LOCALE_PROMISE_VALUE:21,NOT_SUPPORT_LOCALE_ASYNC_FUNCTION:22,NOT_SUPPORT_LOCALE_TYPE:23},Om=24;function Nn(s){return _r(s,null,void 0)}function Ll(s,n){return n.locale!=null?tp(n.locale):tp(s.locale)}let Tt;function tp(s){if(as(s))return s;if(Ls(s)){if(s.resolvedOnce&&Tt!=null)return Tt;if(s.constructor.name==="Function"){const n=s();if(jm(n))throw Nn(In.NOT_SUPPORT_LOCALE_PROMISE_VALUE);return Tt=n}else throw Nn(In.NOT_SUPPORT_LOCALE_ASYNC_FUNCTION)}else throw Nn(In.NOT_SUPPORT_LOCALE_TYPE)}function Fm(s,n,a){return[...new Set([a,...Is(n)?n:Cs(n)?Object.keys(n):as(n)?[n]:[a]])]}function $m(s,n,a){const e=as(a)?a:Ke,t=s;t.__localeChainCache||(t.__localeChainCache=new Map);let l=t.__localeChainCache.get(e);if(!l){l=[];let o=[a];for(;Is(o);)o=lp(l,o,n);const p=Is(n)||!vs(n)?n:n.default?n.default:null;o=as(p)?[p]:p,Is(o)&&lp(l,o,!1),t.__localeChainCache.set(e,l)}return l}function lp(s,n,a){let e=!0;for(let t=0;t<n.length&&Os(e);t++){const l=n[t];as(l)&&(e=Im(s,n[t],a))}return e}function Im(s,n,a){let e;const t=n.split("-");do{const l=t.join("-");e=Nm(s,l,a),t.splice(-1,1)}while(t.length&&e===!0);return e}function Nm(s,n,a){let e=!1;if(!s.includes(n)&&(e=!0,n)){e=n[n.length-1]!=="!";const t=n.replace(/!/g,"");s.push(t),(Is(a)||vs(a))&&a[t]&&(e=a[t])}return e}function Bm(s,n){return Cs(s)?s[n]:null}const zm="12.0.0-alpha.3",it=-1,Ke="en-US",op="",pp=s=>`${s.charAt(0).toLocaleUpperCase()}${s.substr(1)}`;function qm(){return{upper:(s,n)=>n==="text"&&as(s)?s.toUpperCase():n==="vnode"&&Cs(s)&&"__v_isVNode"in s?s.children.toUpperCase():s,lower:(s,n)=>n==="text"&&as(s)?s.toLowerCase():n==="vnode"&&Cs(s)&&"__v_isVNode"in s?s.children.toLowerCase():s,capitalize:(s,n)=>n==="text"&&as(s)?pp(s):n==="vnode"&&Cs(s)&&"__v_isVNode"in s?pp(s.children):s}}let Hm,vr=null;const Vm=s=>{vr=s},Gm=()=>vr;let wr=null;const cp=s=>{wr=s},Um=()=>wr;let rp=0;function Wm(s={}){const n=Ls(s.onWarn)?s.onWarn:bm,a=as(s.version)?s.version:zm,e=as(s.locale)||Ls(s.locale)?s.locale:Ke,t=Ls(e)?Ke:e,l=Is(s.fallbackLocale)||vs(s.fallbackLocale)||as(s.fallbackLocale)||s.fallbackLocale===!1?s.fallbackLocale:t,o=vs(s.messages)?s.messages:Rt(t),p=vs(s.datetimeFormats)?s.datetimeFormats:Rt(t),c=vs(s.numberFormats)?s.numberFormats:Rt(t),u=sn(Ps(),s.modifiers,qm()),r=s.pluralRules||Ps(),i=Ls(s.missing)?s.missing:null,d=Os(s.missingWarn)||We(s.missingWarn)?s.missingWarn:!0,h=Os(s.fallbackWarn)||We(s.fallbackWarn)?s.fallbackWarn:!0,j=!!s.fallbackFormat,m=!!s.unresolving,P=Ls(s.postTranslation)?s.postTranslation:null,C=vs(s.processor)?s.processor:null,b=Os(s.warnHtmlMessage)?s.warnHtmlMessage:!0,x=!!s.escapeParameter,f=Ls(s.messageCompiler)?s.messageCompiler:Hm,k=Ls(s.messageResolver)?s.messageResolver:Bm,w=Ls(s.localeFallbacker)?s.localeFallbacker:Fm,y=Cs(s.fallbackContext)?s.fallbackContext:void 0,M=s,S=Cs(M.__datetimeFormatters)?M.__datetimeFormatters:new Map,F=Cs(M.__numberFormatters)?M.__numberFormatters:new Map,Z=Cs(M.__meta)?M.__meta:{};rp++;const B={version:a,cid:rp,locale:e,fallbackLocale:l,messages:o,modifiers:u,pluralRules:r,missing:i,missingWarn:d,fallbackWarn:h,fallbackFormat:j,unresolving:m,postTranslation:P,processor:C,warnHtmlMessage:b,escapeParameter:x,messageCompiler:f,messageResolver:k,localeFallbacker:w,fallbackContext:y,onWarn:n,__meta:Z};return B.datetimeFormats=p,B.numberFormats=c,B.__datetimeFormatters=S,B.__numberFormatters=F,__INTLIFY_PROD_DEVTOOLS__&&Mm(B,a,Z),B}const Rt=s=>({[s]:Ps()});function Ol(s,n,a,e,t){const{missing:l,onWarn:o}=s;if(l!==null){const p=l(s,a,n,t);return as(p)?p:n}else return n}function Na(s,n,a){const e=s;e.__localeChainCache=new Map,s.localeFallbacker(s,a,n)}function Km(s,n){return s===n?!1:s.split("-")[0]===n.split("-")[0]}function Ym(s,n){const a=n.indexOf(s);if(a===-1)return!1;for(let e=a+1;e<n.length;e++)if(Km(s,n[e]))return!0;return!1}function ip(s,...n){const{datetimeFormats:a,unresolving:e,fallbackLocale:t,onWarn:l,localeFallbacker:o}=s,{__datetimeFormatters:p}=s,[c,u,r,i]=Qt(...n),d=Os(r.missingWarn)?r.missingWarn:s.missingWarn;Os(r.fallbackWarn)?r.fallbackWarn:s.fallbackWarn;const h=!!r.part,j=Ll(s,r),m=o(s,t,j);if(!as(c)||c==="")return new Intl.DateTimeFormat(j,i).format(u);let P={},C,b=null;const x="datetime format";for(let w=0;w<m.length&&(C=m[w],P=a[C]||{},b=P[c],!vs(b));w++)Ol(s,c,C,d,x);if(!vs(b)||!as(C))return e?it:c;let f=`${C}__${c}`;rt(i)||(f=`${f}__${JSON.stringify(i)}`);let k=p.get(f);return k||(k=new Intl.DateTimeFormat(C,sn({},b,i)),p.set(f,k)),h?k.formatToParts(u):k.format(u)}const xr=["localeMatcher","weekday","era","year","month","day","hour","minute","second","timeZoneName","formatMatcher","hour12","timeZone","dateStyle","timeStyle","calendar","dayPeriod","numberingSystem","hourCycle","fractionalSecondDigits"];function Qt(...s){const[n,a,e,t]=s,l=Ps();let o=Ps(),p;if(as(n)){const c=n.match(/(\d{4}-\d{2}-\d{2})(T|\s)?(.*)/);if(!c)throw Nn(In.INVALID_ISO_DATE_ARGUMENT);const u=c[3]?c[3].trim().startsWith("T")?`${c[1].trim()}${c[3].trim()}`:`${c[1].trim()}T${c[3].trim()}`:c[1].trim();p=new Date(u);try{p.toISOString()}catch{throw Nn(In.INVALID_ISO_DATE_ARGUMENT)}}else if(hm(n)){if(isNaN(n.getTime()))throw Nn(In.INVALID_DATE_ARGUMENT);p=n}else if(Zs(n))p=n;else throw Nn(In.INVALID_ARGUMENT);return as(a)?l.key=a:vs(a)&&Object.keys(a).forEach(c=>{xr.includes(c)?o[c]=a[c]:l[c]=a[c]}),as(e)?l.locale=e:vs(e)&&(o=e),vs(t)&&(o=t),[l.key||"",p,l,o]}function up(s,n,a){const e=s;for(const t in a){const l=`${n}__${t}`;e.__datetimeFormatters.has(l)&&e.__datetimeFormatters.delete(l)}}function hp(s,...n){const{numberFormats:a,unresolving:e,fallbackLocale:t,onWarn:l,localeFallbacker:o}=s,{__numberFormatters:p}=s,[c,u,r,i]=Zt(...n),d=Os(r.missingWarn)?r.missingWarn:s.missingWarn;Os(r.fallbackWarn)?r.fallbackWarn:s.fallbackWarn;const h=!!r.part,j=Ll(s,r),m=o(s,t,j);if(!as(c)||c==="")return new Intl.NumberFormat(j,i).format(u);let P={},C,b=null;const x="number format";for(let w=0;w<m.length&&(C=m[w],P=a[C]||{},b=P[c],!vs(b));w++)Ol(s,c,C,d,x);if(!vs(b)||!as(C))return e?it:c;let f=`${C}__${c}`;rt(i)||(f=`${f}__${JSON.stringify(i)}`);let k=p.get(f);return k||(k=new Intl.NumberFormat(C,sn({},b,i)),p.set(f,k)),h?k.formatToParts(u):k.format(u)}const Cr=["localeMatcher","style","currency","currencyDisplay","currencySign","useGrouping","minimumIntegerDigits","minimumFractionDigits","maximumFractionDigits","minimumSignificantDigits","maximumSignificantDigits","compactDisplay","notation","signDisplay","unit","unitDisplay","roundingMode","roundingPriority","roundingIncrement","trailingZeroDisplay"];function Zt(...s){const[n,a,e,t]=s,l=Ps();let o=Ps();if(!Zs(n))throw Nn(In.INVALID_ARGUMENT);const p=n;return as(a)?l.key=a:vs(a)&&Object.keys(a).forEach(c=>{Cr.includes(c)?o[c]=a[c]:l[c]=a[c]}),as(e)?l.locale=e:vs(e)&&(o=e),vs(t)&&(o=t),[l.key||"",p,l,o]}function dp(s,n,a){const e=s;for(const t in a){const l=`${n}__${t}`;e.__numberFormatters.has(l)&&e.__numberFormatters.delete(l)}}const Xm=s=>s,Jm=s=>"",Qm="text",Zm=s=>s.length===0?"":fm(s),sj=gm;function mp(s,n){return s=Math.abs(s),n===2?s?s>1?1:0:1:s?Math.min(s,2):0}function nj(s){const n=Zs(s.pluralIndex)?s.pluralIndex:-1;return s.named&&(Zs(s.named.count)||Zs(s.named.n))?Zs(s.named.count)?s.named.count:Zs(s.named.n)?s.named.n:n:n}function aj(s,n){n.count||(n.count=s),n.n||(n.n=s)}function ej(s={}){const n=s.locale,a=nj(s),e=Cs(s.pluralRules)&&as(n)&&Ls(s.pluralRules[n])?s.pluralRules[n]:mp,t=Cs(s.pluralRules)&&as(n)&&Ls(s.pluralRules[n])?mp:void 0,l=C=>C[e(a,C.length,t)],o=s.list||[],p=C=>o[C],c=s.named||Ps();Zs(s.pluralIndex)&&aj(a,c);const u=C=>c[C];function r(C,b){const x=Ls(s.messages)?s.messages(C,!!b):Cs(s.messages)?s.messages[C]:!1;return x||(s.parent?s.parent.message(C):Jm)}const i=C=>s.modifiers?s.modifiers[C]:Xm,d=vs(s.processor)&&Ls(s.processor.normalize)?s.processor.normalize:Zm,h=vs(s.processor)&&Ls(s.processor.interpolate)?s.processor.interpolate:sj,j=vs(s.processor)&&as(s.processor.type)?s.processor.type:Qm,P={list:p,named:u,plural:l,linked:(C,...b)=>{const[x,f]=b;let k="text",w="";b.length===1?Cs(x)?(w=x.modifier||w,k=x.type||k):as(x)&&(w=x||w):b.length===2&&(as(x)&&(w=x||w),as(f)&&(k=f||k));const y=r(C,!0)(P),M=k==="vnode"&&Is(y)&&w?y[0]:y;return w?i(w)(M,k):M},message:r,type:j,interpolate:h,normalize:d,values:sn(Ps(),o,c)};return P}const jp=()=>"",fn=s=>Ls(s);function gp(s,...n){const{fallbackFormat:a,postTranslation:e,unresolving:t,messageCompiler:l,fallbackLocale:o,messages:p}=s,[c,u]=sl(...n),r=Os(u.missingWarn)?u.missingWarn:s.missingWarn,i=Os(u.fallbackWarn)?u.fallbackWarn:s.fallbackWarn,d=Os(u.escapeParameter)?u.escapeParameter:s.escapeParameter,h=!!u.resolvedMessage,j=as(u.default)||Os(u.default)?Os(u.default)?l?c:()=>c:u.default:a?l?c:()=>c:null,m=a||j!=null&&(as(j)||Ls(j)),P=Ll(s,u);d&&tj(u);let[C,b,x]=h?[c,P,p[P]||Ps()]:kr(s,c,P,o,i,r),f=C,k=c;if(!h&&!(as(f)||Bn(f)||fn(f))&&m&&(f=j,k=f),!h&&(!(as(f)||Bn(f)||fn(f))||!as(b)))return t?it:c;let w=!1;const y=()=>{w=!0},M=fn(f)?f:Sr(s,c,b,f,k,y);if(w)return f;const S=pj(s,b,x,u),F=ej(S),Z=lj(s,M,F),B=e?e(Z,c):Z;if(__INTLIFY_PROD_DEVTOOLS__){const O={timestamp:Date.now(),key:as(c)?c:fn(f)?f.key:"",locale:b||(fn(f)?f.locale:""),format:as(f)?f:fn(f)?f.source:"",message:B};O.meta=sn({},s.__meta,Gm()||{}),Dm(O)}return B}function tj(s){Is(s.list)?s.list=s.list.map(n=>as(n)?ep(n):n):Cs(s.named)&&Object.keys(s.named).forEach(n=>{as(s.named[n])&&(s.named[n]=ep(s.named[n]))})}function kr(s,n,a,e,t,l){const{messages:o,onWarn:p,messageResolver:c,localeFallbacker:u}=s,r=u(s,e,a);let i=Ps(),d,h=null;const j="translate";for(let m=0;m<r.length&&(d=r[m],i=o[d]||Ps(),(h=c(i,n))===null&&(h=i[n]),!(as(h)||Bn(h)||fn(h)));m++)if(!Ym(d,r)){const P=Ol(s,n,d,l,j);P!==n&&(h=P)}return[h,d,i]}function Sr(s,n,a,e,t,l){const{messageCompiler:o,warnHtmlMessage:p}=s;if(fn(e)){const u=e;return u.locale=u.locale||a,u.key=u.key||n,u}if(o==null){const u=(()=>e);return u.locale=a,u.key=n,u}const c=o(e,oj(s,a,t,e,p,l));return c.locale=a,c.key=n,c.source=e,c}function lj(s,n,a){return n(a)}function sl(...s){const[n,a,e]=s,t=Ps();if(!as(n)&&!Zs(n)&&!fn(n)&&!Bn(n))throw Nn(In.INVALID_ARGUMENT);const l=Zs(n)?String(n):(fn(n),n);return Zs(a)?t.plural=a:as(a)?t.default=a:vs(a)&&!rt(a)?t.named=a:Is(a)&&(t.list=a),Zs(e)?t.plural=e:as(e)?t.default=e:vs(e)&&sn(t,e),[l,t]}function oj(s,n,a,e,t,l){return{locale:n,key:a,warnHtmlMessage:t,onError:o=>{throw l&&l(o),o},onCacheKey:o=>im(n,a,o)}}function pj(s,n,a,e){const{modifiers:t,pluralRules:l,messageResolver:o,fallbackLocale:p,fallbackWarn:c,missingWarn:u,fallbackContext:r}=s,d={locale:n,modifiers:t,pluralRules:l,messages:(h,j)=>{let m=o(a,h);if(m==null&&(r||j)){const[,,P]=kr(r||s,h,n,p,c,u);m=o(P,h)}if(as(m)||Bn(m)){let P=!1;const b=Sr(s,h,n,m,h,()=>{P=!0});return P?jp:b}else return fn(m)?m:jp}};return s.processor&&(d.processor=s.processor),e.list&&(d.list=e.list),e.named&&(d.named=e.named),Zs(e.plural)&&(d.pluralIndex=e.plural),d}ym();/*!
  * vue-i18n-core v12.0.0-alpha.3
  * (c) 2016-present kazuya kawaguchi and contributors
  * Released under the MIT License.
  */const cj="12.0.0-alpha.3";function rj(){typeof __VUE_I18N_FULL_INSTALL__!="boolean"&&(Pa().__VUE_I18N_FULL_INSTALL__=!0),typeof __INTLIFY_DROP_MESSAGE_COMPILER__!="boolean"&&(Pa().__INTLIFY_DROP_MESSAGE_COMPILER__=!1),typeof __INTLIFY_PROD_DEVTOOLS__!="boolean"&&(Pa().__INTLIFY_PROD_DEVTOOLS__=!1)}const Hn={UNEXPECTED_RETURN_TYPE:Om,INVALID_ARGUMENT:25,MUST_BE_CALL_SETUP_TOP:26,NOT_INSTALLED:27,NOT_INSTALLED_WITH_PROVIDE:31,UNEXPECTED_ERROR:32};function ea(s,...n){return _r(s,null,void 0)}const nl=pa("__translateVNode"),al=pa("__datetimeParts"),el=pa("__numberParts"),ij=pa("__setPluralRules"),uj=pa("__injectWithOption"),tl=pa("__dispose");function ce(s){if(!Cs(s)||Bn(s))return s;for(const n in s)if(fa(s,n))if(!n.includes("."))Cs(s[n])&&ce(s[n]);else{const a=n.split("."),e=a.length-1;let t=s,l=!1;for(let o=0;o<e;o++){if(a[o]==="__proto__")throw new Error(`unsafe key: ${a[o]}`);if(a[o]in t||(t[a[o]]=Ps()),!Cs(t[a[o]])){l=!0;break}t=t[a[o]]}if(l||(Bn(t)?Rm.includes(a[e])||delete s[n]:(t[a[e]]=s[n],delete s[n])),!Bn(t)){const o=t[a[e]];Cs(o)&&ce(o)}}return s}function Pr(s,n){const{messages:a,__i18n:e,messageResolver:t,flatJson:l}=n,o=vs(a)?a:Is(e)?Ps():{[s]:Ps()};if(Is(e)&&e.forEach(p=>{if("locale"in p&&"resource"in p){const{locale:c,resource:u}=p;c?(o[c]=o[c]||Ps(),Le(u,o[c])):Le(u,o)}else as(p)&&Le(JSON.parse(p),o)}),t==null&&l)for(const p in o)fa(o,p)&&ce(o[p]);return o}function Ar(s){return s.type}function hj(s,n,a){let e=Cs(n.messages)?n.messages:Ps();"__i18nGlobal"in a&&(e=Pr(s.locale.value,{messages:e,__i18n:a.__i18nGlobal}));const t=Object.keys(e);t.length&&t.forEach(l=>{s.mergeLocaleMessage(l,e[l])});{if(Cs(n.datetimeFormats)){const l=Object.keys(n.datetimeFormats);l.length&&l.forEach(o=>{s.mergeDateTimeFormat(o,n.datetimeFormats[o])})}if(Cs(n.numberFormats)){const l=Object.keys(n.numberFormats);l.length&&l.forEach(o=>{s.mergeNumberFormat(o,n.numberFormats[o])})}}}function fp(s){return bs(me,null,s,0)}const bp="__INTLIFY_META__",_p=()=>[],dj=()=>!1;let yp=0;function vp(s){return((n,a,e,t)=>s(a,e,Gn()||void 0,t))}const mj=()=>{const s=Gn();let n=null;return s&&(n=Ar(s)[bp])?{[bp]:n}:null};function Tr(s={}){const{__root:n,__injectWithOption:a}=s,e=n===void 0,t=s.flatJson,l=np?ns:_l;let o=Os(s.inheritLocale)?s.inheritLocale:!0;const p=l(n&&o?n.locale.value:as(s.locale)?s.locale:Ke),c=l(n&&o?n.fallbackLocale.value:as(s.fallbackLocale)||Is(s.fallbackLocale)||vs(s.fallbackLocale)||s.fallbackLocale===!1?s.fallbackLocale:p.value),u=l(Pr(p.value,s)),r=l(vs(s.datetimeFormats)?s.datetimeFormats:{[p.value]:{}}),i=l(vs(s.numberFormats)?s.numberFormats:{[p.value]:{}});let d=n?n.missingWarn:Os(s.missingWarn)||We(s.missingWarn)?s.missingWarn:!0,h=n?n.fallbackWarn:Os(s.fallbackWarn)||We(s.fallbackWarn)?s.fallbackWarn:!0,j=n?n.fallbackRoot:Os(s.fallbackRoot)?s.fallbackRoot:!0,m=!!s.fallbackFormat,P=Ls(s.missing)?s.missing:null,C=Ls(s.missing)?vp(s.missing):null,b=Ls(s.postTranslation)?s.postTranslation:null,x=n?n.warnHtmlMessage:Os(s.warnHtmlMessage)?s.warnHtmlMessage:!0,f=!!s.escapeParameter;const k=n?n.modifiers:vs(s.modifiers)?s.modifiers:{};let w=s.pluralRules||n&&n.pluralRules,y;y=(()=>{e&&cp(null);const T={version:cj,locale:p.value,fallbackLocale:c.value,messages:u.value,modifiers:k,pluralRules:w,missing:C===null?void 0:C,missingWarn:d,fallbackWarn:h,fallbackFormat:m,unresolving:!0,postTranslation:b===null?void 0:b,warnHtmlMessage:x,escapeParameter:f,messageResolver:s.messageResolver,messageCompiler:s.messageCompiler,__meta:{framework:"vue"}};T.datetimeFormats=r.value,T.numberFormats=i.value,T.__datetimeFormatters=vs(y)?y.__datetimeFormatters:void 0,T.__numberFormatters=vs(y)?y.__numberFormatters:void 0;const D=Wm(T);return e&&cp(D),D})(),Na(y,p.value,c.value);function S(){return[p.value,c.value,u.value,r.value,i.value]}const F=ss({get:()=>p.value,set:T=>{y.locale=T,p.value=T}}),Z=ss({get:()=>c.value,set:T=>{y.fallbackLocale=T,c.value=T,Na(y,p.value,T)}}),B=ss(()=>u.value),O=ss(()=>Object.keys(u.value).sort()),V=ss(()=>r.value),us=ss(()=>i.value);function os(){return Ls(b)?b:null}function Y(T){b=T,y.postTranslation=T}function hs(){return P}function Fs(T){T!==null&&(C=vp(T)),P=T,y.missing=C}const $s=(T,D,ts,cs,Ds,qs)=>{S();let Us;try{__INTLIFY_PROD_DEVTOOLS__,e||(y.fallbackContext=n?Um():void 0),Us=T(y)}finally{__INTLIFY_PROD_DEVTOOLS__,e||(y.fallbackContext=void 0)}if(ts!=="translate exists"&&Zs(Us)&&Us===it||ts==="translate exists"&&!Us){const[jn,ba]=D();return n&&j?cs(n):Ds(jn)}else{if(qs(Us))return Us;throw ea(Hn.UNEXPECTED_RETURN_TYPE)}};function Hs(...T){return $s(D=>Reflect.apply(gp,null,[D,...T]),()=>sl(...T),"translate",D=>Reflect.apply(D.t,D,[...T]),D=>D,D=>as(D))}function Vs(...T){const[D,ts,cs]=T;if(cs&&!Cs(cs))throw ea(Hn.INVALID_ARGUMENT);return Hs(D,ts,sn({resolvedMessage:!0},cs||{}))}function es(...T){return $s(D=>Reflect.apply(ip,null,[D,...T]),()=>Qt(...T),"datetime format",D=>Reflect.apply(D.d,D,[...T]),()=>op,D=>as(D)||Is(D))}function ms(...T){return $s(D=>Reflect.apply(hp,null,[D,...T]),()=>Zt(...T),"number format",D=>Reflect.apply(D.n,D,[...T]),()=>op,D=>as(D)||Is(D))}function _s(T){return T.map(D=>as(D)||Zs(D)||Os(D)?fp(String(D)):D)}const R={normalize:_s,interpolate:T=>T,type:"vnode"};function H(...T){return $s(D=>{let ts;const cs=D;try{cs.processor=R,ts=Reflect.apply(gp,null,[cs,...T])}finally{cs.processor=null}return ts},()=>sl(...T),"translate",D=>D[nl](...T),D=>[fp(D)],D=>Is(D))}function q(...T){return $s(D=>Reflect.apply(hp,null,[D,...T]),()=>Zt(...T),"number format",D=>D[el](...T),_p,D=>as(D)||Is(D))}function X(...T){return $s(D=>Reflect.apply(ip,null,[D,...T]),()=>Qt(...T),"datetime format",D=>D[al](...T),_p,D=>as(D)||Is(D))}function gs(T){w=T,y.pluralRules=w}function g(T,D){return $s(()=>{if(!T)return!1;const ts=as(D)?D:p.value,cs=$(ts),Ds=y.messageResolver(cs,T);return Bn(Ds)||fn(Ds)||as(Ds)},()=>[T],"translate exists",ts=>Reflect.apply(ts.te,ts,[T,D]),dj,ts=>Os(ts))}function _(T){let D=null;const ts=$m(y,c.value,p.value);for(let cs=0;cs<ts.length;cs++){const Ds=u.value[ts[cs]]||{},qs=y.messageResolver(Ds,T);if(qs!=null){D=qs;break}}return D}function A(T){const D=_(T);return D??(n?n.tm(T)||{}:{})}function $(T){return u.value[T]||{}}function N(T,D){if(t){const ts={[T]:D};for(const cs in ts)fa(ts,cs)&&ce(ts[cs]);D=ts[T]}u.value[T]=D,y.messages=u.value}function I(T,D){u.value[T]=u.value[T]||{};const ts={[T]:D};if(t)for(const cs in ts)fa(ts,cs)&&ce(ts[cs]);D=ts[T],Le(D,u.value[T]),y.messages=u.value}function K(T){return r.value[T]||{}}function W(T,D){r.value[T]=D,y.datetimeFormats=r.value,up(y,T,D)}function U(T,D){r.value[T]=sn(r.value[T]||{},D),y.datetimeFormats=r.value,up(y,T,D)}function z(T){return i.value[T]||{}}function ps(T,D){i.value[T]=D,y.numberFormats=i.value,dp(y,T,D)}function J(T,D){i.value[T]=sn(i.value[T]||{},D),y.numberFormats=i.value,dp(y,T,D)}yp++,n&&np&&(Bs(n.locale,T=>{o&&(p.value=T,y.locale=T,Na(y,p.value,c.value))}),Bs(n.fallbackLocale,T=>{o&&(c.value=T,y.fallbackLocale=T,Na(y,p.value,c.value))}));const Q={id:yp,locale:F,fallbackLocale:Z,get inheritLocale(){return o},set inheritLocale(T){o=T,T&&n&&(p.value=n.locale.value,c.value=n.fallbackLocale.value,Na(y,p.value,c.value))},availableLocales:O,messages:B,get modifiers(){return k},get pluralRules(){return w||{}},get isGlobal(){return e},get missingWarn(){return d},set missingWarn(T){d=T,y.missingWarn=d},get fallbackWarn(){return h},set fallbackWarn(T){h=T,y.fallbackWarn=h},get fallbackRoot(){return j},set fallbackRoot(T){j=T},get fallbackFormat(){return m},set fallbackFormat(T){m=T,y.fallbackFormat=m},get warnHtmlMessage(){return x},set warnHtmlMessage(T){x=T,y.warnHtmlMessage=T},get escapeParameter(){return f},set escapeParameter(T){f=T,y.escapeParameter=T},t:Hs,getLocaleMessage:$,setLocaleMessage:N,mergeLocaleMessage:I,getPostTranslationHandler:os,setPostTranslationHandler:Y,getMissingHandler:hs,setMissingHandler:Fs,[ij]:gs};return Q.datetimeFormats=V,Q.numberFormats=us,Q.rt=Vs,Q.te=g,Q.tm=A,Q.d=es,Q.n=ms,Q.getDateTimeFormat=K,Q.setDateTimeFormat=W,Q.mergeDateTimeFormat=U,Q.getNumberFormat=z,Q.setNumberFormat=ps,Q.mergeNumberFormat=J,Q[uj]=a,Q[nl]=H,Q[al]=X,Q[el]=q,Q}const Fl={tag:{type:[String,Object]},locale:{type:String},scope:{type:String,validator:s=>s==="parent"||s==="global",default:"parent"},i18n:{type:Object}};function jj({slots:s},n){return n.length===1&&n[0]==="default"?(s.default?s.default():[]).reduce((e,t)=>[...e,...t.type===js?t.children:[t]],[]):n.reduce((a,e)=>{const t=s[e];return t&&(a[e]=t()),a},Ps())}function Rr(){return js}function gj(s){return Is(s)&&!as(s[0])}function Er(s,n,a,e){const{slots:t,attrs:l}=n;return()=>{const o={part:!0};let p=Ps();s.locale&&(o.locale=s.locale),as(s.format)?o.key=s.format:Cs(s.format)&&(as(s.format.key)&&(o.key=s.format.key),p=Object.keys(s.format).reduce((d,h)=>a.includes(h)?sn(Ps(),d,{[h]:s.format[h]}):d,Ps()));const c=e(s.value,o,p);let u=[o.key];Is(c)?u=c.map((d,h)=>{const j=t[d.type],m=j?j({[d.type]:d.value,index:h,parts:c}):[d.value];return gj(m)&&(m[0].key=`${d.type}-${h}`),m}):as(c)&&(u=[c]);const r=sn(Ps(),l),i=as(s.tag)||Cs(s.tag)?s.tag:Rr();return ge(i,r,u)}}const fj=zs({name:"i18n-d",props:sn({value:{type:[Number,Date],required:!0},format:{type:[String,Object]}},Fl),setup(s,n){const a=s.i18n||Xs({useScope:s.scope,__useComponent:!0});return Er(s,n,xr,(...e)=>a[al](...e))}}),wp=fj,bj=zs({name:"i18n-n",props:sn({value:{type:Number,required:!0},format:{type:[String,Object]}},Fl),setup(s,n){const a=s.i18n||Xs({useScope:s.scope,__useComponent:!0});return Er(s,n,Cr,(...e)=>a[el](...e))}}),xp=bj,_j=zs({name:"i18n-t",props:sn({},{keypath:{type:String,required:!0},plural:{type:[Number,String],validator:s=>Zs(s)||!isNaN(s)}},Fl),setup(s,n){const{slots:a,attrs:e}=n,t=s.i18n||Xs({useScope:s.scope,__useComponent:!0});return()=>{const l=Object.keys(a).filter(i=>i[0]!=="_"),o=Ps();s.locale&&(o.locale=s.locale),s.plural!==void 0&&(o.plural=as(s.plural)?+s.plural:s.plural);const p=jj(n,l),c=t[nl](s.keypath,p,o),u=sn(Ps(),e),r=as(s.tag)||Cs(s.tag)?s.tag:Rr();return ge(r,u,c)}}}),Cp=_j;function yj(s,...n){const a=vs(n[0])?n[0]:{};(Os(a.globalInstall)?a.globalInstall:!0)&&([Cp.name,"I18nT"].forEach(t=>s.component(t,Cp)),[xp.name,"I18nN"].forEach(t=>s.component(t,xp)),[wp.name,"I18nD"].forEach(t=>s.component(t,wp)))}const vj=pa("global-vue-i18n");function wj(s={}){const n=Os(s.globalInjection)?s.globalInjection:!0,a=new Map,[e,t]=xj(s),l=pa("");function o(r){return a.get(r)||null}function p(r,i){a.set(r,i)}function c(r){a.delete(r)}const u={async install(r,...i){if(r.__VUE_I18N_SYMBOL__=l,r.provide(r.__VUE_I18N_SYMBOL__,u),vs(i[0])){const j=i[0];u.__composerExtend=j.__composerExtend}let d=null;n&&(d=Ej(r,u.global)),__VUE_I18N_FULL_INSTALL__&&yj(r,...i);const h=r.unmount;r.unmount=()=>{d&&d(),u.dispose(),h()}},get global(){return t},dispose(){e.stop()},__instances:a,__getInstance:o,__setInstance:p,__deleteInstance:c};return u}function Xs(s={}){const n=Gn();if(n==null)throw ea(Hn.MUST_BE_CALL_SETUP_TOP);if(!n.isCE&&n.appContext.app!=null&&!n.appContext.app.__VUE_I18N_SYMBOL__)throw ea(Hn.NOT_INSTALLED);const a=Cj(n),e=Sj(a),t=Ar(n),l=kj(s,t);if(l==="global")return hj(e,s,t),e;if(l==="parent"){let c=Pj(a,n,s.__useComponent);return c==null&&(c=e),c}const o=a;let p=o.__getInstance(n);if(p==null){const c=sn({},s);"__i18n"in t&&(c.__i18n=t.__i18n),e&&(c.__root=e),p=Tr(c),o.__composerExtend&&(p[tl]=o.__composerExtend(p)),Tj(o,n,p),o.__setInstance(n,p)}return p}function xj(s){const n=il(),a=n.run(()=>Tr(s));if(a==null)throw ea(Hn.UNEXPECTED_ERROR);return[n,a]}function Cj(s){const n=on(s.isCE?vj:s.appContext.app.__VUE_I18N_SYMBOL__);if(!n)throw ea(s.isCE?Hn.NOT_INSTALLED_WITH_PROVIDE:Hn.UNEXPECTED_ERROR);return n}function kj(s,n){return rt(s)?"__i18n"in n?"local":"global":s.useScope?s.useScope:"local"}function Sj(s){return s.global}function Pj(s,n,a=!1){let e=null;const t=n.root;let l=Aj(n,a);for(;l!=null&&(e=s.__getInstance(l),!(e!=null||t===l));)l=l.parent;return e}function Aj(s,n=!1){return s==null?null:n&&s.vnode.ctx||s.parent}function Tj(s,n,a){dn(()=>{},n),wl(()=>{const e=a;s.__deleteInstance(n);const t=e[tl];t&&(t(),delete e[tl])},n)}const Rj=["locale","fallbackLocale","availableLocales"],kp=["t","rt","d","n","tm","te"];function Ej(s,n){const a=Object.create(null);return Rj.forEach(t=>{const l=Object.getOwnPropertyDescriptor(n,t);if(!l)throw ea(Hn.UNEXPECTED_ERROR);const o=Es(l.value)?{get(){return l.value.value},set(p){l.value.value=p}}:{get(){return l.get&&l.get()}};Object.defineProperty(a,t,o)}),s.config.globalProperties.$i18n=a,kp.forEach(t=>{const l=Object.getOwnPropertyDescriptor(n,t);if(!l||!l.value)throw ea(Hn.UNEXPECTED_ERROR);Object.defineProperty(s.config.globalProperties,`$${t}`,l)}),()=>{delete s.config.globalProperties.$i18n,kp.forEach(t=>{delete s.config.globalProperties[`$${t}`]})}}rj();if(__INTLIFY_PROD_DEVTOOLS__){const s=Pa();s.__INTLIFY__=!0,Em(s.__INTLIFY_DEVTOOLS_GLOBAL_HOOK__)}const Sp=s=>{if(typeof window>"u"||typeof document>"u")return;const n=document.documentElement;if(s==="auto"){const a=window.matchMedia("(prefers-color-scheme: dark)").matches;n.setAttribute("data-bs-theme",a?"dark":"light"),localStorage.removeItem("theme")}else n.setAttribute("data-bs-theme",s),localStorage.setItem("theme",s)},Mj={tags:"标签",articles:"文章",words:"字数",prevPage:"上一页",nextPage:"下一页",pagination:"分页导航",search:"搜索",theme:"主题",language:"语言",menu:"菜单",close:"关闭",searchPlaceholder:"搜索文章标题或正文...",searchNoResults:"未找到匹配的文章",searchUnavailable:"搜索索引加载失败",searchEnterHint:"Enter 打开",searchEscHint:"Esc 关闭",searchResultsLabel:"条结果",home:"首页",categories:"分类",resources:"资源",about:"关于",greeting:"你好，我是",greetingPrefix:"//",avatar:"头像",developer:"开发者",wordUnit:"字",recentPosts:"最近文章",notes:"笔记",projects:"项目",topics:"课题",seeMore:"查看更多",countPosts:"{count} 篇",countProjects:"{count} 个项目",countTopics:"{count} 个课题",countWords:"{count} 字",copyCode:"复制代码",copyFailed:"复制失败",copyArticle:"复制文章",copyTable:"复制表格",copied:"已复制",anchorHeading:"置顶当前标题",uncategorized:"未分类",updatedAt:"更新于",postReadingTime:"{minutes} 分钟",articleReadingTime:"阅读约 {minutes} 分钟",tableOfContents:"本页目录",openToc:"打开目录",backToTop:"回到顶部",resourceSubtitle:"生物信息学与结构生物学领域常用工具",resourceItems:"项",experience:"经历",introduction:"自我介绍",contact:"联系我",thanks:"感谢您的关注！",designedByPrefix:"由",designedBySuffix:"设计",clearFilter:"清除筛选",backToArticle:"返回文章"},Dj={tags:"Tags",articles:"Articles",words:"Words",prevPage:"Previous",nextPage:"Next",pagination:"Pagination",search:"Search",theme:"Theme",language:"Language",menu:"Menu",close:"Close",searchPlaceholder:"Search articles by title or content...",searchNoResults:"No matching articles found",searchUnavailable:"Search index failed to load",searchEnterHint:"Enter to open",searchEscHint:"Esc to close",searchResultsLabel:"results",home:"Home",categories:"Categories",resources:"Resources",about:"About",greeting:"Hello, I'm",greetingPrefix:"//",avatar:"Avatar",developer:"Developer",wordUnit:"words",recentPosts:"Recent Posts",notes:"Notes",projects:"Projects",topics:"Topics",seeMore:"See More",countPosts:"{count} posts",countProjects:"{count} projects",countTopics:"{count} topics",countWords:"{count} words",copyCode:"Copy code",copyFailed:"Copy failed",copyArticle:"Copy article",copyTable:"Copy table",copied:"Copied",anchorHeading:"Anchor to heading",uncategorized:"Uncategorized",updatedAt:"Updated at",postReadingTime:"{minutes} min",articleReadingTime:"Reading about {minutes} minutes",tableOfContents:"Table of Contents",openToc:"Open table of contents",backToTop:"Back to top",resourceSubtitle:"Commonly used tools in bioinformatics and structural biology",resourceItems:"items",experience:"Experience",introduction:"Introduction",contact:"Contact Me",thanks:"Thank you for your attention!",designedByPrefix:"Designed by",designedBySuffix:"",clearFilter:"Clear filter",backToArticle:"Back to article"},en={url:"https://zorrooz.github.io",author:"zorrooz",email:"zorrooz@163.com",github:"https://github.com/zorrooz",startYear:2025},_n=["zh-CN","en-US"],Qn=["/zh","/en"],Et=["auto","light","dark"],Mr=s=>s==="en"?_n[1]:_n[0],Dr=typeof window<"u"&&localStorage.getItem("locale")||_n[0];typeof document<"u"&&(document.documentElement.lang=Dr);const Ye=wj({locale:Dr,fallbackLocale:_n[0],messages:{[_n[0]]:Mj,[_n[1]]:Dj},globalInjection:!0}),$l=rm("app",()=>{const s=ns(typeof window<"u"&&localStorage.getItem("theme")||"auto"),n=ns(typeof window<"u"&&localStorage.getItem("locale")||_n[0]);return{theme:s,locale:n,setLocale:o=>{n.value=o,Ye.global.locale.value=o,typeof window<"u"&&localStorage.setItem("locale",o),typeof document<"u"&&(document.documentElement.lang=o)},toggleTheme:()=>{const p=(Et.indexOf(s.value)+1)%Et.length;s.value=Et[p],Sp(s.value)},initTheme:()=>{Sp(s.value)},initLocale:()=>{if(typeof window<"u"){const o=window.location.pathname.match(/^\/(zh|en)(\/|$)/);o&&(n.value=Mr(o[1]))}Ye.global.locale.value=n.value,typeof document<"u"&&(document.documentElement.lang=n.value)}}}),$n=s=>`${Ye.global.locale.value===_n[1]?Qn[1]:Qn[0]}${s.startsWith("/")?s:`/${s}`}`,Lr=`Hello, I am zorrooz, currently focused on structural analysis and functional studies of biological macromolecules, while also exploring computational biology.
`,Or=[{year:"2021 – 2025",title:"Lanzhou University · B.S. in Bioinformatics",desc:"Basic programming training and foundational knowledge in biology"},{year:"2025 – present",title:"Lanzhou University · M.S. in Biology",desc:"Structural and functional studies of biological macromolecules"}],Fr=[{title:"Common Tools",items:[{name:"Python · R · Bash · Git",desc:"Primary languages and version control tools for daily research and development"},{name:"Nextflow · Snakemake",desc:"Bioinformatics pipeline orchestration and workflow management"},{name:"RELION · cryoSPARC",desc:"Cryo-EM single-particle reconstruction and 3D structure validation"},{name:"VS Code · Ubuntu · WSL",desc:"Development environment and terminal workflow"}]},{title:"Fields",items:[{name:"Structural Biology",desc:"Cryo-EM"},{name:"Computational Biology",desc:"Binder design and virtual screening"},{name:"Programming",desc:"R · Python · C · Perl · Bash · JavaScript"}]}],$r=[{label:"Email",value:"zorrooz@163.com",link:"mailto:zorrooz@163.com",icon:"fas fa-envelope"},{label:"GitHub",value:"zorrooz",link:"https://github.com/zorrooz",icon:"fab fa-github"}],Lj={introduction:Lr,experience:Or,section:Fr,contacts:$r},Oj=Object.freeze(Object.defineProperty({__proto__:null,contacts:$r,default:Lj,experience:Or,introduction:Lr,section:Fr},Symbol.toStringTag,{value:"Module"})),Ir=`你好，我是 zorrooz，当前专注于生物大分子的结构解析与功能研究，同时涉猎计算生物学。
`,Nr=[{year:"2021 – 2025",title:"兰州大学 本科 · 生物信息学",desc:"基本编程训练与生物学基础知识"},{year:"2025 – 至今",title:"兰州大学 硕士 · 生物学",desc:"生物大分子的结构与功能研究"}],Br=[{title:"常用工具",items:[{name:"Python · R · Bash · Git",desc:"日常研究与开发的主力语言与版本控制工具"},{name:"Nextflow · Snakemake",desc:"生信流水线编排与工作流管理"},{name:"RELION · cryoSPARC",desc:"冷冻电镜单颗粒重构与三维结构验证"},{name:"VS Code · Ubuntu · WSL",desc:"开发环境与终端工作流"}]},{title:"领域",items:[{name:"结构生物学",desc:"冷冻电镜"},{name:"计算生物学",desc:"binder 设计与虚拟筛选"},{name:"编程",desc:"R · Python · C · Perl · Bash · JavaScript"}]}],zr=[{label:"邮箱",value:"zorrooz@163.com",link:"mailto:zorrooz@163.com",icon:"fas fa-envelope"},{label:"GitHub",value:"zorrooz",link:"https://github.com/zorrooz",icon:"fab fa-github"}],Fj={introduction:Ir,experience:Nr,section:Br,contacts:zr},$j=Object.freeze(Object.defineProperty({__proto__:null,contacts:zr,default:Fj,experience:Nr,introduction:Ir,section:Br},Symbol.toStringTag,{value:"Module"})),Ij=JSON.parse('[{"title":"Notes","items":[{"name":"ComputationalBiology","title":"Computational Biology","desc":"Mainstream tools for protein design and virtual screening","tags":["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology","Virtual Screening","Molecular Docking","AutoDock Vina"],"stats":{"postsCount":2,"totalWords":1211,"latestDate":"2026-08-04"},"categories":[{"key":"protein-design","title":"Protein Design","articles":[{"title":"Mainstream Tools for Protein Design","articleUrl":"/article/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en","wordCount":618,"tags":["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology"]}],"stats":{"postsCount":1,"totalWords":618,"latestDate":"2026-08-04"}},{"key":"virtual-screening","title":"Virtual Screening","articles":[{"title":"Mainstream Tools for Virtual Screening","articleUrl":"/article/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools-en","wordCount":593,"tags":["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"]}],"stats":{"postsCount":1,"totalWords":593,"latestDate":"2026-08-04"}}],"root":"/article/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en"},{"name":"Programming","title":"Programming","desc":"Detailed tutorials on R, Python, Linux, and Bash programming","tags":["Bash","Shell","Scripting","Tutorial","Linux","Command Line","Python","Advanced","Getting Started","NumPy","Pandas","Data Processing","R Language","ggplot2","Data Visualization","tidyverse","dplyr"],"stats":{"postsCount":8,"totalWords":1972,"latestDate":"2026-08-04"},"categories":[{"key":"python","title":"Python","articles":[{"title":"Python Advanced: Functions, Classes, and Modules","articleUrl":"/article/notes/Programming/python/python-advanced/python-advanced-en","wordCount":244,"tags":["Python","Advanced","Tutorial"]},{"title":"Introduction to Python Programming: Environment, Syntax, and Data Types","articleUrl":"/article/notes/Programming/python/python-basics/python-basics-en","wordCount":333,"tags":["Python","Getting Started","Tutorial"]},{"title":"Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas","articleUrl":"/article/notes/Programming/python/python-data/python-data-en","wordCount":314,"tags":["Python","NumPy","Pandas","Data Processing","Tutorial"]}],"stats":{"postsCount":3,"totalWords":891,"latestDate":"2026-08-04"}},{"key":"r","title":"R","articles":[{"title":"R Language Primer: Data Structures, Vectorization, and Functions","articleUrl":"/article/notes/Programming/r/r-basics/r-basics-en","wordCount":256,"tags":["R Language","Getting Started","Tutorial"]},{"title":"R ggplot2 Data Visualization: Grammar of Graphics and Common Charts","articleUrl":"/article/notes/Programming/r/r-ggplot2/r-ggplot2-en","wordCount":157,"tags":["R Language","ggplot2","Data Visualization","Tutorial"]},{"title":"R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning","articleUrl":"/article/notes/Programming/r/r-tidyverse/r-tidyverse-en","wordCount":245,"tags":["R Language","tidyverse","dplyr","Tutorial"]}],"stats":{"postsCount":3,"totalWords":658,"latestDate":"2026-08-04"}},{"key":"linux","title":"Linux","articles":[{"title":"Linux Command Line Basics: File System, Permissions, and Text Processing","articleUrl":"/article/notes/Programming/linux/linux-basics/linux-basics-en","wordCount":238,"tags":["Linux","Command Line","Tutorial"]}],"stats":{"postsCount":1,"totalWords":238,"latestDate":"2026-08-04"}},{"key":"bash","title":"Bash","articles":[{"title":"Bash Programming: Variables, Conditionals, Loops, and Practical Scripts","articleUrl":"/article/notes/Programming/bash/bash-scripting/bash-scripting-en","wordCount":185,"tags":["Bash","Shell","Scripting","Tutorial"]}],"stats":{"postsCount":1,"totalWords":185,"latestDate":"2026-08-04"}}],"root":"/article/notes/Programming/python/python-advanced/python-advanced-en"},{"name":"StructuralBiology","title":"Structural Biology","desc":"Cryo-EM data processing, visualization of biomacromolecular structures, and atomic modeling","tags":["Cryo-Electron Microscopy","cryo-EM","Review","Data Processing","RELION","Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial","Structure Visualization","PyMOL","ChimeraX"],"stats":{"postsCount":4,"totalWords":2670,"latestDate":"2026-08-04"},"categories":[{"key":"cryoem","title":"Cryo-EM","articles":[{"title":"Review of Cryo-EM Technology","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en","wordCount":977,"tags":["Cryo-Electron Microscopy","cryo-EM","Review"]},{"title":"Cryo-EM Single Particle Analysis: Full Data Processing Workflow","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow-en","wordCount":765,"tags":["Cryo-Electron Microscopy","cryo-EM","Data Processing","RELION"]}],"stats":{"postsCount":2,"totalWords":1742,"latestDate":"2026-08-04"}},{"key":"visualization","title":"Structural Visualization","articles":[{"title":"Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice","articleUrl":"/article/notes/StructuralBiology/visualization/structure-visualization/structure-visualization-en","wordCount":362,"tags":["Structure Visualization","PyMOL","ChimeraX","Tutorial"]}],"stats":{"postsCount":1,"totalWords":362,"latestDate":"2026-08-04"}},{"key":"modeling","title":"Atomic Modeling","articles":[{"title":"Atomic Modeling and Refinement","articleUrl":"/article/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling-en","wordCount":566,"tags":["Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial"]}],"stats":{"postsCount":1,"totalWords":566,"latestDate":"2026-08-04"}}],"root":"/article/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en"}]},{"title":"Projects","items":[{"name":"Plotit","title":"plotit","desc":"Declarative plotting package for R — a verb-prefix API based on ggplot2, in early development stage","tags":["R","ggplot2","Data Visualization","Declarative API"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-12-30"},"categories":[],"root":"","github":"https://github.com/zorrooz/plotit","status":"active","language":"R","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/plotit"},{"name":"ZorroozBlog","title":"zorrooz.github.io","desc":"Personal blog system — built with Vue 3 + Vite, supporting bilingual Chinese/English, theme switching, and Markdown article management","tags":["Vue","Vite","Blog","Markdown"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-08-29"},"categories":[],"root":"","github":"https://github.com/zorrooz/zorrooz.github.io","status":"active","language":"Vue","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/"}]},{"title":"Topics","items":[{"name":"StructureOfProteinDemo","title":"Protein Structure Determination Demo Project","desc":"Protein structure placeholder demo","tags":["Protein","Structure","Cryo-EM"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2026-04-12"},"categories":[],"root":"","doi":"10.1234/demo-structure.2026","status":"completed","journal":"Journal of Demo Science","year":2026,"authors":["Demo Author A","Demo Author B"],"url":"https://www.demo-structure.org"}]}]'),Nj=Object.freeze(Object.defineProperty({__proto__:null,default:Ij},Symbol.toStringTag,{value:"Module"})),Bj=JSON.parse('[{"title":"笔记","items":[{"name":"ComputationalBiology","title":"计算生物学","desc":"蛋白质设计与虚拟筛选的主流工具","tags":["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学","虚拟筛选","分子对接","AutoDock Vina"],"stats":{"postsCount":2,"totalWords":1899,"latestDate":"2026-08-04"},"categories":[{"key":"protein-design","title":"蛋白质设计","articles":[{"title":"蛋白质设计主流工具","articleUrl":"/article/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools","wordCount":1002,"tags":["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学"]}],"stats":{"postsCount":1,"totalWords":1002,"latestDate":"2026-08-04"}},{"key":"virtual-screening","title":"虚拟筛选","articles":[{"title":"虚拟筛选主流工具","articleUrl":"/article/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools","wordCount":897,"tags":["虚拟筛选","分子对接","AutoDock Vina","计算生物学"]}],"stats":{"postsCount":1,"totalWords":897,"latestDate":"2026-08-04"}}],"root":"/article/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools"},{"name":"Programming","title":"编程","desc":"R、Python、Linux 与 Bash 编程的详细教程","tags":["Bash","Shell","脚本编程","教程","Linux","命令行","Python","进阶","入门","NumPy","Pandas","数据处理","R语言","ggplot2","数据可视化","tidyverse","dplyr"],"stats":{"postsCount":8,"totalWords":2903,"latestDate":"2026-08-04"},"categories":[{"key":"python","title":"Python","articles":[{"title":"Python 进阶：函数、类与模块","articleUrl":"/article/notes/Programming/python/python-advanced/python-advanced","wordCount":376,"tags":["Python","进阶","教程"]},{"title":"Python 编程入门：环境、语法与数据类型","articleUrl":"/article/notes/Programming/python/python-basics/python-basics","wordCount":495,"tags":["Python","入门","教程"]},{"title":"Python 数据处理实战：文件 IO、正则与 NumPy/Pandas","articleUrl":"/article/notes/Programming/python/python-data/python-data","wordCount":473,"tags":["Python","NumPy","Pandas","数据处理","教程"]}],"stats":{"postsCount":3,"totalWords":1344,"latestDate":"2026-08-04"}},{"key":"r","title":"R","articles":[{"title":"R 语言入门：数据结构、向量化与函数","articleUrl":"/article/notes/Programming/r/r-basics/r-basics","wordCount":403,"tags":["R语言","入门","教程"]},{"title":"R ggplot2 数据可视化：图层语法与常用图表","articleUrl":"/article/notes/Programming/r/r-ggplot2/r-ggplot2","wordCount":245,"tags":["R语言","ggplot2","数据可视化","教程"]},{"title":"R tidyverse 数据操作：dplyr 管道与数据清洗","articleUrl":"/article/notes/Programming/r/r-tidyverse/r-tidyverse","wordCount":324,"tags":["R语言","tidyverse","dplyr","教程"]}],"stats":{"postsCount":3,"totalWords":972,"latestDate":"2026-08-04"}},{"key":"linux","title":"Linux","articles":[{"title":"Linux 命令行基础：文件系统、权限与文本处理","articleUrl":"/article/notes/Programming/linux/linux-basics/linux-basics","wordCount":328,"tags":["Linux","命令行","教程"]}],"stats":{"postsCount":1,"totalWords":328,"latestDate":"2026-08-04"}},{"key":"bash","title":"Bash","articles":[{"title":"Bash 编程：变量、条件、循环与实用脚本","articleUrl":"/article/notes/Programming/bash/bash-scripting/bash-scripting","wordCount":259,"tags":["Bash","Shell","脚本编程","教程"]}],"stats":{"postsCount":1,"totalWords":259,"latestDate":"2026-08-04"}}],"root":"/article/notes/Programming/python/python-advanced/python-advanced"},{"name":"StructuralBiology","title":"结构生物学","desc":"冷冻电镜数据处理、生物大分子结构可视化与原子建模","tags":["冷冻电镜","cryo-EM","综述","数据处理","RELION","原子建模","Coot","phenix","结构精修","教程","结构可视化","PyMOL","ChimeraX"],"stats":{"postsCount":4,"totalWords":3960,"latestDate":"2026-08-04"},"categories":[{"key":"cryoem","title":"冷冻电镜","articles":[{"title":"冷冻电镜技术综述","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview","wordCount":1485,"tags":["冷冻电镜","cryo-EM","综述"]},{"title":"冷冻电镜单颗粒分析：数据处理全流程","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow","wordCount":1111,"tags":["冷冻电镜","cryo-EM","数据处理","RELION"]}],"stats":{"postsCount":2,"totalWords":2596,"latestDate":"2026-08-04"}},{"key":"visualization","title":"结构可视化","articles":[{"title":"生物大分子结构可视化：PyMOL 与 ChimeraX 实战","articleUrl":"/article/notes/StructuralBiology/visualization/structure-visualization/structure-visualization","wordCount":529,"tags":["结构可视化","PyMOL","ChimeraX","教程"]}],"stats":{"postsCount":1,"totalWords":529,"latestDate":"2026-08-04"}},{"key":"modeling","title":"原子建模","articles":[{"title":"原子建模与精修","articleUrl":"/article/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling","wordCount":835,"tags":["原子建模","Coot","phenix","结构精修","教程"]}],"stats":{"postsCount":1,"totalWords":835,"latestDate":"2026-08-04"}}],"root":"/article/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview"}]},{"title":"项目","items":[{"name":"Plotit","title":"plotit","desc":"R 声明式绘图包——基于 ggplot2 的动词前缀 API，早期开发阶段","tags":["R","ggplot2","数据可视化","声明式API"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-12-30"},"categories":[],"root":"","github":"https://github.com/zorrooz/plotit","status":"active","language":"R","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/plotit"},{"name":"ZorroozBlog","title":"zorrooz.github.io","desc":"个人博客系统——基于 Vue 3 + Vite 构建，支持中英双语、主题切换与 Markdown 文章管理","tags":["Vue","Vite","博客","Markdown"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-08-29"},"categories":[],"root":"","github":"https://github.com/zorrooz/zorrooz.github.io","status":"active","language":"Vue","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/"}]},{"title":"课题","items":[{"name":"StructureOfProteinDemo","title":"蛋白质结构解析演示课题","desc":"蛋白质结构占位demo","tags":["Protein","Structure","Cryo-EM"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2026-04-12"},"categories":[],"root":"","doi":"10.1234/demo-structure.2026","status":"completed","journal":"Journal of Demo Science","year":2026,"authors":["Demo Author A","Demo Author B"],"url":"https://www.demo-structure.org"}]}]'),zj=Object.freeze(Object.defineProperty({__proto__:null,default:Bj},Symbol.toStringTag,{value:"Module"})),qj=[{title:"Mainstream Tools for Protein Design",date:"2026-08-04",tags:["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology"],draft:!1,description:"A panorama of protein design tools: Rosetta energy optimization, AlphaFold structure prediction, RFdiffusion generative design, ESMFold fast prediction, and binder design workflow",relativePath:"ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en",wordCount:618,tagCount:5},{title:"Mainstream Tools for Virtual Screening",date:"2026-08-04",tags:["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"],draft:!1,description:"A comprehensive overview of virtual screening tools: docking software such as AutoDock Vina, Glide, and DOCK, receptor/ligand preparation, scoring functions, high-throughput screening workflows, and validation strategies",relativePath:"ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools-en",wordCount:593,tagCount:4},{title:"Bash Programming: Variables, Conditionals, Loops, and Practical Scripts",date:"2026-08-04",tags:["Bash","Shell","Scripting","Tutorial"],draft:!1,description:"A complete Bash scripting tutorial: variables and parameters, conditional statements, loops, functions, arrays, and string processing, with practical scripts for batch processing and bioinformatics",relativePath:"Programming/bash/bash-scripting/bash-scripting-en",wordCount:185,tagCount:4},{title:"Linux Command Line Basics: File System, Permissions, and Text Processing",date:"2026-08-04",tags:["Linux","Command Line","Tutorial"],draft:!1,description:"Core skills for the Linux command line: file system navigation, file operations, permission management, the three musketeers of text processing (grep/sed/awk), and process management",relativePath:"Programming/linux/linux-basics/linux-basics-en",wordCount:238,tagCount:3},{title:"Python Advanced: Functions, Classes, and Modules",date:"2026-08-04",tags:["Python","Advanced","Tutorial"],draft:!1,description:"Deep dive into Python functional and object-oriented programming: parameter passing, lambda, decorators, classes and inheritance, exception handling, modules and packages",relativePath:"Programming/python/python-advanced/python-advanced-en",wordCount:244,tagCount:3},{title:"Introduction to Python Programming: Environment, Syntax, and Data Types",date:"2026-08-04",tags:["Python","Getting Started","Tutorial"],draft:!1,description:"Learn Python from scratch: environment setup, basic syntax, built-in data types, and control flow, with runnable example code.",relativePath:"Programming/python/python-basics/python-basics-en",wordCount:333,tagCount:3},{title:"Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas",date:"2026-08-04",tags:["Python","NumPy","Pandas","Data Processing","Tutorial"],draft:!1,description:"A complete practical tutorial on file I/O, regular expressions, and the scientific computing stack (NumPy/Pandas), covering common data scenarios in bioinformatics",relativePath:"Programming/python/python-data/python-data-en",wordCount:314,tagCount:5},{title:"R Language Primer: Data Structures, Vectorization, and Functions",date:"2026-08-04",tags:["R Language","Getting Started","Tutorial"],draft:!1,description:"Core R fundamentals: environment setup, vectors/matrices/data frames/lists, vectorized operations, control flow, and function definitions",relativePath:"Programming/r/r-basics/r-basics-en",wordCount:256,tagCount:3},{title:"R ggplot2 Data Visualization: Grammar of Graphics and Common Charts",date:"2026-08-04",tags:["R Language","ggplot2","Data Visualization","Tutorial"],draft:!1,description:"Introduction to the ggplot2 Grammar of Graphics: scatter plots, boxplots, histograms, bar charts, and theme customization for publication-ready figures",relativePath:"Programming/r/r-ggplot2/r-ggplot2-en",wordCount:157,tagCount:4},{title:"R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning",date:"2026-08-04",tags:["R Language","tidyverse","dplyr","Tutorial"],draft:!1,description:"Chain dplyr verbs (filter/select/mutate/group_by/summarise/join) with pipes, and use tidyr to complete a full data cleaning workflow.",relativePath:"Programming/r/r-tidyverse/r-tidyverse-en",wordCount:245,tagCount:4},{title:"Review of Cryo-EM Technology",date:"2026-08-04",tags:["Cryo-Electron Microscopy","cryo-EM","Review"],draft:!1,description:"The technological revolution of single-particle cryo-EM: direct electron detectors, the stability revolution, and the advent of the AI era, along with the future landscape of structural biology",relativePath:"StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en",wordCount:977,tagCount:3},{title:"Cryo-EM Single Particle Analysis: Full Data Processing Workflow",date:"2026-08-04",tags:["Cryo-Electron Microscopy","cryo-EM","Data Processing","RELION"],draft:!1,description:"From micrographs to atomic models: the complete workflow of cryo-EM single particle analysis (SPA), covering preprocessing, particle picking, 2D/3D classification, refinement, and resolution assessment",relativePath:"StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow-en",wordCount:765,tagCount:4},{title:"Atomic Modeling and Refinement",date:"2026-08-04",tags:["Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial"],draft:!1,description:"Complete workflow for atomic model building and refinement from cryo-EM density maps: automatic modeling, Coot manual correction, phenix refinement, model quality assessment metrics",relativePath:"StructuralBiology/modeling/atomic-modeling/atomic-modeling-en",wordCount:566,tagCount:5},{title:"Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice",date:"2026-08-04",tags:["Structure Visualization","PyMOL","ChimeraX","Tutorial"],draft:!1,description:"Visualizing protein/nucleic acid structures with PyMOL and UCSF ChimeraX: PDB data retrieval, rendering modes, script batch processing, and cryo-EM density map display",relativePath:"StructuralBiology/visualization/structure-visualization/structure-visualization-en",wordCount:362,tagCount:4}],Hj=Object.freeze(Object.defineProperty({__proto__:null,default:qj},Symbol.toStringTag,{value:"Module"})),Vj=[{title:"蛋白质设计主流工具",date:"2026-08-04",tags:["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学"],draft:!1,description:"蛋白质设计工具全景：Rosetta 能量优化、AlphaFold 结构预测、RFdiffusion 生成式设计、ESMFold 快速预测与 binder 设计工作流",relativePath:"ComputationalBiology/protein-design/protein-design-tools/protein-design-tools",wordCount:1002,tagCount:5},{title:"虚拟筛选主流工具",date:"2026-08-04",tags:["虚拟筛选","分子对接","AutoDock Vina","计算生物学"],draft:!1,description:"虚拟筛选工具全景：AutoDock Vina、Glide、DOCK 等对接软件，受体/配体准备、打分函数、高通量筛选流程与验证策略",relativePath:"ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools",wordCount:897,tagCount:4},{title:"Bash 编程：变量、条件、循环与实用脚本",date:"2026-08-04",tags:["Bash","Shell","脚本编程","教程"],draft:!1,description:"Bash 脚本编程完整教程：变量与参数、条件判断、循环、函数、数组与字符串处理，附批处理与生信实用脚本",relativePath:"Programming/bash/bash-scripting/bash-scripting",wordCount:259,tagCount:4},{title:"Linux 命令行基础：文件系统、权限与文本处理",date:"2026-08-04",tags:["Linux","命令行","教程"],draft:!1,description:"Linux 命令行核心技能：文件系统导航、文件操作、权限管理、文本处理三剑客（grep/sed/awk）与进程管理",relativePath:"Programming/linux/linux-basics/linux-basics",wordCount:328,tagCount:3},{title:"Python 进阶：函数、类与模块",date:"2026-08-04",tags:["Python","进阶","教程"],draft:!1,description:"深入 Python 函数式与面向对象编程：参数传递、lambda、装饰器、类与继承、异常处理、模块与包",relativePath:"Programming/python/python-advanced/python-advanced",wordCount:376,tagCount:3},{title:"Python 编程入门：环境、语法与数据类型",date:"2026-08-04",tags:["Python","入门","教程"],draft:!1,description:"从零开始学习 Python：环境搭建、基础语法、内置数据类型与流程控制，附带可运行的示例代码",relativePath:"Programming/python/python-basics/python-basics",wordCount:495,tagCount:3},{title:"Python 数据处理实战：文件 IO、正则与 NumPy/Pandas",date:"2026-08-04",tags:["Python","NumPy","Pandas","数据处理","教程"],draft:!1,description:"文件读写、正则表达式与科学计算三件套（NumPy/Pandas）的完整实战教程，覆盖生物信息学常见数据场景",relativePath:"Programming/python/python-data/python-data",wordCount:473,tagCount:5},{title:"R 语言入门：数据结构、向量化与函数",date:"2026-08-04",tags:["R语言","入门","教程"],draft:!1,description:"R 语言核心基础：环境配置、向量/矩阵/数据框/列表、向量化运算、控制流与函数定义",relativePath:"Programming/r/r-basics/r-basics",wordCount:403,tagCount:3},{title:"R ggplot2 数据可视化：图层语法与常用图表",date:"2026-08-04",tags:["R语言","ggplot2","数据可视化","教程"],draft:!1,description:"ggplot2 图层语法入门：散点图、箱线图、直方图、柱状图与主题美化，输出出版级图表",relativePath:"Programming/r/r-ggplot2/r-ggplot2",wordCount:245,tagCount:4},{title:"R tidyverse 数据操作：dplyr 管道与数据清洗",date:"2026-08-04",tags:["R语言","tidyverse","dplyr","教程"],draft:!1,description:"用管道符串联 dplyr 动词（filter/select/mutate/group_by/summarise/join），配合 tidyr 完成完整的数据清洗流程",relativePath:"Programming/r/r-tidyverse/r-tidyverse",wordCount:324,tagCount:4},{title:"冷冻电镜技术综述",date:"2026-08-04",tags:["冷冻电镜","cryo-EM","综述"],draft:!1,description:"冷冻电镜单颗粒分析的技术革命：直接电子探测器、稳定性革命与 AI 时代的到来，以及结构生物学的未来图景",relativePath:"StructuralBiology/cryoem/cryoem-overview/cryoem-overview",wordCount:1485,tagCount:3},{title:"冷冻电镜单颗粒分析：数据处理全流程",date:"2026-08-04",tags:["冷冻电镜","cryo-EM","数据处理","RELION"],draft:!1,description:"从显微照片到原子模型：冷冻电镜单颗粒分析（SPA）完整流程，涵盖预处理、颗粒挑选、2D/3D 分类、精修与分辨率评估",relativePath:"StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow",wordCount:1111,tagCount:4},{title:"原子建模与精修",date:"2026-08-04",tags:["原子建模","Coot","phenix","结构精修","教程"],draft:!1,description:"冷冻电镜密度图的原子模型搭建与精修全流程：自动建模、Coot 手动修正、phenix 精修、模型质量评估指标",relativePath:"StructuralBiology/modeling/atomic-modeling/atomic-modeling",wordCount:835,tagCount:5},{title:"生物大分子结构可视化：PyMOL 与 ChimeraX 实战",date:"2026-08-04",tags:["结构可视化","PyMOL","ChimeraX","教程"],draft:!1,description:"用 PyMOL 与 UCSF ChimeraX 可视化蛋白质/核酸结构：PDB 数据获取、渲染模式、脚本批处理与冷冻电镜密度图展示",relativePath:"StructuralBiology/visualization/structure-visualization/structure-visualization",wordCount:529,tagCount:4}],Gj=Object.freeze(Object.defineProperty({__proto__:null,default:Vj},Symbol.toStringTag,{value:"Module"})),Uj=[{id:1,no:1,title:"Mainstream Tools for Protein Design",date:"2026-08-04",category:["ComputationalBiology","protein-design"],tags:["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology"],preview:"A panorama of protein design tools: Rosetta energy optimization, AlphaFold structure prediction, RFdiffusion generative design, ESMFold fast prediction, and binder design workflow",wordCount:618},{id:2,no:2,title:"Mainstream Tools for Virtual Screening",date:"2026-08-04",category:["ComputationalBiology","virtual-screening"],tags:["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"],preview:"A comprehensive overview of virtual screening tools: docking software such as AutoDock Vina, Glide, and DOCK, receptor/ligand preparation, scoring functions, high-throughput screening workflows, and validation strategies",wordCount:593},{id:3,no:3,title:"Bash Programming: Variables, Conditionals, Loops, and Practical Scripts",date:"2026-08-04",category:["Programming","bash"],tags:["Bash","Shell","Scripting","Tutorial"],preview:"A complete Bash scripting tutorial: variables and parameters, conditional statements, loops, functions, arrays, and string processing, with practical scripts for batch processing and bioinformatics",wordCount:185},{id:4,no:4,title:"Linux Command Line Basics: File System, Permissions, and Text Processing",date:"2026-08-04",category:["Programming","linux"],tags:["Linux","Command Line","Tutorial"],preview:"Core skills for the Linux command line: file system navigation, file operations, permission management, the three musketeers of text processing (grep/sed/awk), and process management",wordCount:238},{id:5,no:5,title:"Python Advanced: Functions, Classes, and Modules",date:"2026-08-04",category:["Programming","python"],tags:["Python","Advanced","Tutorial"],preview:"Deep dive into Python functional and object-oriented programming: parameter passing, lambda, decorators, classes and inheritance, exception handling, modules and packages",wordCount:244},{id:6,no:6,title:"Introduction to Python Programming: Environment, Syntax, and Data Types",date:"2026-08-04",category:["Programming","python"],tags:["Python","Getting Started","Tutorial"],preview:"Learn Python from scratch: environment setup, basic syntax, built-in data types, and control flow, with runnable example code.",wordCount:333},{id:7,no:7,title:"Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas",date:"2026-08-04",category:["Programming","python"],tags:["Python","NumPy","Pandas","Data Processing","Tutorial"],preview:"A complete practical tutorial on file I/O, regular expressions, and the scientific computing stack (NumPy/Pandas), covering common data scenarios in bioinformatics",wordCount:314},{id:8,no:8,title:"R Language Primer: Data Structures, Vectorization, and Functions",date:"2026-08-04",category:["Programming","r"],tags:["R Language","Getting Started","Tutorial"],preview:"Core R fundamentals: environment setup, vectors/matrices/data frames/lists, vectorized operations, control flow, and function definitions",wordCount:256},{id:9,no:9,title:"R ggplot2 Data Visualization: Grammar of Graphics and Common Charts",date:"2026-08-04",category:["Programming","r"],tags:["R Language","ggplot2","Data Visualization","Tutorial"],preview:"Introduction to the ggplot2 Grammar of Graphics: scatter plots, boxplots, histograms, bar charts, and theme customization for publication-ready figures",wordCount:157},{id:10,no:10,title:"R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning",date:"2026-08-04",category:["Programming","r"],tags:["R Language","tidyverse","dplyr","Tutorial"],preview:"Chain dplyr verbs (filter/select/mutate/group_by/summarise/join) with pipes, and use tidyr to complete a full data cleaning workflow.",wordCount:245},{id:11,no:11,title:"Review of Cryo-EM Technology",date:"2026-08-04",category:["StructuralBiology","cryoem"],tags:["Cryo-Electron Microscopy","cryo-EM","Review"],preview:"The technological revolution of single-particle cryo-EM: direct electron detectors, the stability revolution, and the advent of the AI era, along with the future landscape of structural biology",wordCount:977},{id:12,no:12,title:"Cryo-EM Single Particle Analysis: Full Data Processing Workflow",date:"2026-08-04",category:["StructuralBiology","cryoem"],tags:["Cryo-Electron Microscopy","cryo-EM","Data Processing","RELION"],preview:"From micrographs to atomic models: the complete workflow of cryo-EM single particle analysis (SPA), covering preprocessing, particle picking, 2D/3D classification, refinement, and resolution assessment",wordCount:765},{id:13,no:13,title:"Atomic Modeling and Refinement",date:"2026-08-04",category:["StructuralBiology","modeling"],tags:["Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial"],preview:"Complete workflow for atomic model building and refinement from cryo-EM density maps: automatic modeling, Coot manual correction, phenix refinement, model quality assessment metrics",wordCount:566},{id:14,no:14,title:"Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice",date:"2026-08-04",category:["StructuralBiology","visualization"],tags:["Structure Visualization","PyMOL","ChimeraX","Tutorial"],preview:"Visualizing protein/nucleic acid structures with PyMOL and UCSF ChimeraX: PDB data retrieval, rendering modes, script batch processing, and cryo-EM density map display",wordCount:362}],Wj=Object.freeze(Object.defineProperty({__proto__:null,default:Uj},Symbol.toStringTag,{value:"Module"})),Kj=[{id:1,no:1,title:"蛋白质设计主流工具",date:"2026-08-04",category:["计算生物学","protein-design"],tags:["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学"],preview:"蛋白质设计工具全景：Rosetta 能量优化、AlphaFold 结构预测、RFdiffusion 生成式设计、ESMFold 快速预测与 binder 设计工作流",wordCount:1002},{id:2,no:2,title:"虚拟筛选主流工具",date:"2026-08-04",category:["计算生物学","virtual-screening"],tags:["虚拟筛选","分子对接","AutoDock Vina","计算生物学"],preview:"虚拟筛选工具全景：AutoDock Vina、Glide、DOCK 等对接软件，受体/配体准备、打分函数、高通量筛选流程与验证策略",wordCount:897},{id:3,no:3,title:"Bash 编程：变量、条件、循环与实用脚本",date:"2026-08-04",category:["编程","bash"],tags:["Bash","Shell","脚本编程","教程"],preview:"Bash 脚本编程完整教程：变量与参数、条件判断、循环、函数、数组与字符串处理，附批处理与生信实用脚本",wordCount:259},{id:4,no:4,title:"Linux 命令行基础：文件系统、权限与文本处理",date:"2026-08-04",category:["编程","linux"],tags:["Linux","命令行","教程"],preview:"Linux 命令行核心技能：文件系统导航、文件操作、权限管理、文本处理三剑客（grep/sed/awk）与进程管理",wordCount:328},{id:5,no:5,title:"Python 进阶：函数、类与模块",date:"2026-08-04",category:["编程","python"],tags:["Python","进阶","教程"],preview:"深入 Python 函数式与面向对象编程：参数传递、lambda、装饰器、类与继承、异常处理、模块与包",wordCount:376},{id:6,no:6,title:"Python 编程入门：环境、语法与数据类型",date:"2026-08-04",category:["编程","python"],tags:["Python","入门","教程"],preview:"从零开始学习 Python：环境搭建、基础语法、内置数据类型与流程控制，附带可运行的示例代码",wordCount:495},{id:7,no:7,title:"Python 数据处理实战：文件 IO、正则与 NumPy/Pandas",date:"2026-08-04",category:["编程","python"],tags:["Python","NumPy","Pandas","数据处理","教程"],preview:"文件读写、正则表达式与科学计算三件套（NumPy/Pandas）的完整实战教程，覆盖生物信息学常见数据场景",wordCount:473},{id:8,no:8,title:"R 语言入门：数据结构、向量化与函数",date:"2026-08-04",category:["编程","r"],tags:["R语言","入门","教程"],preview:"R 语言核心基础：环境配置、向量/矩阵/数据框/列表、向量化运算、控制流与函数定义",wordCount:403},{id:9,no:9,title:"R ggplot2 数据可视化：图层语法与常用图表",date:"2026-08-04",category:["编程","r"],tags:["R语言","ggplot2","数据可视化","教程"],preview:"ggplot2 图层语法入门：散点图、箱线图、直方图、柱状图与主题美化，输出出版级图表",wordCount:245},{id:10,no:10,title:"R tidyverse 数据操作：dplyr 管道与数据清洗",date:"2026-08-04",category:["编程","r"],tags:["R语言","tidyverse","dplyr","教程"],preview:"用管道符串联 dplyr 动词（filter/select/mutate/group_by/summarise/join），配合 tidyr 完成完整的数据清洗流程",wordCount:324},{id:11,no:11,title:"冷冻电镜技术综述",date:"2026-08-04",category:["结构生物学","cryoem"],tags:["冷冻电镜","cryo-EM","综述"],preview:"冷冻电镜单颗粒分析的技术革命：直接电子探测器、稳定性革命与 AI 时代的到来，以及结构生物学的未来图景",wordCount:1485},{id:12,no:12,title:"冷冻电镜单颗粒分析：数据处理全流程",date:"2026-08-04",category:["结构生物学","cryoem"],tags:["冷冻电镜","cryo-EM","数据处理","RELION"],preview:"从显微照片到原子模型：冷冻电镜单颗粒分析（SPA）完整流程，涵盖预处理、颗粒挑选、2D/3D 分类、精修与分辨率评估",wordCount:1111},{id:13,no:13,title:"原子建模与精修",date:"2026-08-04",category:["结构生物学","modeling"],tags:["原子建模","Coot","phenix","结构精修","教程"],preview:"冷冻电镜密度图的原子模型搭建与精修全流程：自动建模、Coot 手动修正、phenix 精修、模型质量评估指标",wordCount:835},{id:14,no:14,title:"生物大分子结构可视化：PyMOL 与 ChimeraX 实战",date:"2026-08-04",category:["结构生物学","visualization"],tags:["结构可视化","PyMOL","ChimeraX","教程"],preview:"用 PyMOL 与 UCSF ChimeraX 可视化蛋白质/核酸结构：PDB 数据获取、渲染模式、脚本批处理与冷冻电镜密度图展示",wordCount:529}],Yj=Object.freeze(Object.defineProperty({__proto__:null,default:Kj},Symbol.toStringTag,{value:"Module"})),Xj=[{name:"Plotit",title:"plotit",desc:"Declarative plotting package for R — a verb-prefix API based on ggplot2, in early development stage",date:"2025-12-30",tags:["R","ggplot2","Data Visualization","Declarative API"],github:"https://github.com/zorrooz/plotit",url:"https://zorrooz.github.io/plotit",status:"active",language:"R",stars:0,license:"MIT",version:""},{name:"ZorroozBlog",title:"zorrooz.github.io",desc:"Personal blog system — built with Vue 3 + Vite, supporting bilingual Chinese/English, theme switching, and Markdown article management",date:"2025-08-29",tags:["Vue","Vite","Blog","Markdown"],github:"https://github.com/zorrooz/zorrooz.github.io",url:"https://zorrooz.github.io/",status:"active",language:"Vue",stars:0,license:"MIT",version:""}],Jj=Object.freeze(Object.defineProperty({__proto__:null,default:Xj},Symbol.toStringTag,{value:"Module"})),Qj=[{name:"Plotit",title:"plotit",desc:"R 声明式绘图包——基于 ggplot2 的动词前缀 API，早期开发阶段",date:"2025-12-30",tags:["R","ggplot2","数据可视化","声明式API"],github:"https://github.com/zorrooz/plotit",url:"https://zorrooz.github.io/plotit",status:"active",language:"R",stars:0,license:"MIT",version:""},{name:"ZorroozBlog",title:"zorrooz.github.io",desc:"个人博客系统——基于 Vue 3 + Vite 构建，支持中英双语、主题切换与 Markdown 文章管理",date:"2025-08-29",tags:["Vue","Vite","博客","Markdown"],github:"https://github.com/zorrooz/zorrooz.github.io",url:"https://zorrooz.github.io/",status:"active",language:"Vue",stars:0,license:"MIT",version:""}],Zj=Object.freeze(Object.defineProperty({__proto__:null,default:Qj},Symbol.toStringTag,{value:"Module"})),sg=JSON.parse(`[{"title":"Data Analysis","children":[{"title":"Data Exploration","children":[{"title":"R Language","items":[{"name":"tidyverse","url":"https://www.tidyverse.org/","desc":"A suite of R packages for data science, covering data import, cleaning, transformation, and visualization"},{"name":"data.table","url":"https://rdatatable.gitlab.io/data.table/","desc":"High-performance data manipulation package, extremely fast for processing data frames with millions of rows"},{"name":"readr","url":"https://readr.tidyverse.org/","desc":"A tool for fast reading of CSV/TSV and other tabular files, with automatic column type inference"}]},{"title":"Python","items":[{"name":"NumPy","url":"https://numpy.org/","desc":"The fundamental library for scientific computing in Python, providing multidimensional arrays and vectorized operations"},{"name":"Pandas","url":"https://pandas.pydata.org/","desc":"A powerful data analysis library based on Python, providing efficient data manipulation and processing capabilities"},{"name":"Polars","url":"https://www.pola.rs/","desc":"A high-performance data processing library based on Apache Arrow, supporting multiple language interfaces"},{"name":"Jupyter","url":"https://jupyter.org/","desc":"An interactive notebook environment that supports mixing code, text, and visualizations"}]}]},{"title":"Data Visualization","children":[{"title":"R Language","items":[{"name":"ggplot2","url":"https://ggplot2.tidyverse.org/","desc":"A powerful data visualization package in R, based on the grammar of graphics"},{"name":"plotly (R)","url":"https://plotly.com/r/","desc":"The R interface to the interactive chart library, capable of generating web-interactive graphics"},{"name":"patchwork","url":"https://patchwork.data-imaginist.com/","desc":"A multi-panel composition tool that easily combines ggplot graphics using + and / symbols"},{"name":"RColorBrewer","url":"https://cran.r-project.org/web/packages/RColorBrewer/","desc":"Publication-quality color schemes, providing carefully designed discrete palettes"}]},{"title":"Python","items":[{"name":"Matplotlib","url":"https://matplotlib.org/","desc":"A widely used plotting library in Python, supporting many chart types"},{"name":"Seaborn","url":"https://seaborn.pydata.org/","desc":"A statistical visualization library based on Matplotlib, with built-in attractive themes and statistical plots"},{"name":"Plotly","url":"https://plotly.com/python/","desc":"An interactive visualization library supporting scatter, 3D, maps, and many other chart types"},{"name":"Altair","url":"https://altair-viz.github.io/","desc":"A declarative statistical visualization library based on Vega-Lite"}]}]},{"title":"Statistical Analysis","children":[{"title":"R Language","items":[{"name":"R Stats","url":"https://www.r-project.org/","desc":"The built-in statistical analysis functionality in R, covering a wide range of statistical methods and models"},{"name":"rstatix","url":"https://rpkgs.datanovia.com/rstatix/","desc":"Pipe-friendly wrappers for statistical tests (t-tests, ANOVA, rank-sum tests, etc.)"},{"name":"broom","url":"https://broom.tidymodels.org/","desc":"A tool for tidying statistical model outputs into clean data frames"},{"name":"lme4","url":"https://cran.r-project.org/web/packages/lme4/","desc":"Linear and nonlinear mixed-effects models, suitable for repeated-measures data"}]},{"title":"Python","items":[{"name":"SciPy Stats","url":"https://scipy.org/","desc":"The statistics module in the Python scientific computing library SciPy, providing many statistical distributions and test methods"},{"name":"statsmodels","url":"https://www.statsmodels.org/","desc":"A statistical modeling library supporting regression, time series, hypothesis testing, and more"},{"name":"pingouin","url":"https://pingouin-stats.org/","desc":"A simple and user-friendly statistical testing library covering common parametric and nonparametric tests"}]}]},{"title":"Machine Learning","children":[{"title":"Python","items":[{"name":"scikit-learn","url":"https://scikit-learn.org/","desc":"A machine learning library based on Python, providing a rich set of algorithms and tools"},{"name":"PyTorch","url":"https://pytorch.org/","desc":"A deep learning framework supporting dynamic computation graphs and efficient tensor operations"},{"name":"TensorFlow","url":"https://www.tensorflow.org/","desc":"A deep learning framework with a mature ecosystem, supporting production-level deployment"},{"name":"XGBoost","url":"https://xgboost.readthedocs.io/","desc":"A gradient boosting tree library, the top choice for tabular data competitions and engineering"},{"name":"LightGBM","url":"https://lightgbm.readthedocs.io/","desc":"A gradient boosting framework from Microsoft, with fast training speed and low memory usage"}]},{"title":"R Language","items":[{"name":"tidymodels","url":"https://www.tidymodels.org/","desc":"A unified modeling framework in R, covering data preprocessing, modeling, and evaluation"}]}]}]},{"title":"Omics","children":[{"title":"Data Platforms","children":[{"title":"Sequence Databases","items":[{"name":"NCBI","url":"https://www.ncbi.nlm.nih.gov/","desc":"The National Center for Biotechnology Information in the United States, providing important databases such as GenBank"},{"name":"EMBL-EBI","url":"https://www.ebi.ac.uk/","desc":"The European Molecular Biology Laboratory, providing a variety of bioinformatics resources and tools"},{"name":"DDBJ","url":"https://www.ddbj.nig.ac.jp/","desc":"The DNA Data Bank of Japan, providing storage and access for nucleic acid and protein sequence data"},{"name":"UniProt","url":"https://www.uniprot.org/","desc":"The most comprehensive database for protein sequences and functional annotation"},{"name":"Ensembl","url":"https://www.ensembl.org/","desc":"A genome annotation and comparative genomics database for vertebrates"},{"name":"UCSC Genome Browser","url":"https://genome.ucsc.edu/","desc":"A genome browser supporting visualization of multiple species genomes and custom tracks"}]},{"title":"Sequencing Data","items":[{"name":"SRA","url":"https://www.ncbi.nlm.nih.gov/sra","desc":"The NCBI Sequence Read Archive, storing raw high-throughput sequencing data"},{"name":"GEO","url":"https://www.ncbi.nlm.nih.gov/geo/","desc":"A gene expression database containing microarray and sequencing expression data"},{"name":"ENA","url":"https://www.ebi.ac.uk/ena","desc":"The European Nucleotide Archive, Europe's storage center for raw sequencing data"}]},{"title":"Protein Interactions and Pathways","items":[{"name":"STRING","url":"https://string-db.org/","desc":"A protein-protein interaction database supporting multiple species"},{"name":"KEGG","url":"https://www.kegg.jp/","desc":"A comprehensive database of pathways, genes, and compounds"},{"name":"GO","url":"http://geneontology.org/","desc":"Gene Ontology, providing a standardized system for gene function annotation"}]}]},{"title":"Analysis Software","children":[{"title":"Command Line","items":[{"name":"FastQC","url":"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/","desc":"A quality assessment tool for sequencing data, generating visual QC reports"},{"name":"Trimmomatic","url":"http://www.usadellab.org/cms/?page=trimmomatic","desc":"A read trimming tool for sequencing data, removing adapters and low-quality bases"},{"name":"cutadapt","url":"https://cutadapt.readthedocs.io/","desc":"A command-line adapter removal tool supporting multiple sequencing platforms"},{"name":"BWA","url":"https://github.com/lh3/bwa","desc":"A short-read alignment tool, the standard choice for genome alignment"},{"name":"STAR","url":"https://github.com/alexdobin/STAR","desc":"An ultra-fast splice-aware aligner for RNA-seq"},{"name":"HISAT2","url":"https://daehwankimlab.github.io/hisat2/","desc":"A graph-based index alignment tool for RNA-seq"},{"name":"SAMtools","url":"https://www.htslib.org/","desc":"A toolset for processing SAM/BAM files, the standard for post-alignment processing"},{"name":"bcftools","url":"https://samtools.github.io/bcftools/","desc":"A tool for processing VCF/BCF variant files, supporting filtering and statistics"},{"name":"GATK","url":"https://gatk.broadinstitute.org/","desc":"The Genome Analysis Toolkit, the industry standard for variant detection"},{"name":"MultiQC","url":"https://multiqc.info/","desc":"Aggregates multiple QC reports into a single interactive HTML report"}]},{"title":"Python","items":[{"name":"Biopython","url":"https://biopython.org/","desc":"The core Python library for bioinformatics, covering sequences, structures, and file parsing"},{"name":"scikit-bio","url":"https://scikit-bio.org/","desc":"A Python library for bioinformatics, providing various tools for biological data processing and analysis"},{"name":"scanpy","url":"https://scanpy.readthedocs.io/","desc":"A Python library for single-cell transcriptomics analysis, deeply integrated with AnnData"}]},{"title":"R Language","items":[{"name":"Bioconductor","url":"https://bioconductor.org/","desc":"A collection of bioinformatics software packages based on R, covering genomics, transcriptomics, and other fields"},{"name":"Seurat","url":"https://satijalab.org/seurat/","desc":"A mainstream R package for single-cell transcriptomics analysis, providing a complete analysis workflow"},{"name":"DESeq2","url":"https://bioconductor.org/packages/release/bioc/html/DESeq2.html","desc":"The standard tool for differential expression analysis based on the negative binomial model"},{"name":"edgeR","url":"https://bioconductor.org/packages/release/bioc/html/edgeR.html","desc":"A differential expression analysis tool based on empirical Bayes"},{"name":"limma","url":"https://bioconductor.org/packages/release/bioc/html/limma.html","desc":"A differential analysis tool using linear models, applicable to both microarray and sequencing data"}]}]}]},{"title":"Cryo-EM Structure Determination","children":[{"title":"Software Tools","children":[{"title":"3D Reconstruction","items":[{"name":"RELION","url":"https://www3.mrc-lmb.cam.ac.uk/groups/scheres/index.php","desc":"A Bayesian approach-based software for cryo-EM 3D reconstruction, widely used for high-resolution structure determination"},{"name":"cryoSPARC","url":"https://cryoem.slac.stanford.edu/","desc":"A high-performance cryo-EM data processing platform with an intuitive user interface"},{"name":"cisTEM","url":"https://cistem.org/","desc":"A free, open-source cryo-EM processing software supporting a complete single-particle workflow"},{"name":"SPHIRE","url":"https://sphire.mpg.de/","desc":"A modular cryo-EM processing suite supporting high-resolution reconstruction"}]},{"title":"Preprocessing","items":[{"name":"MotionCor2","url":"https://emcore.ucsf.edu/ucsf-software","desc":"A frame-by-frame motion correction tool that eliminates sample drift caused by the electron beam"},{"name":"CTFFIND4","url":"https://grigoriefflab.umassmed.edu/ctffind4","desc":"A classic tool for CTF estimation, evaluating defocus and astigmatism"},{"name":"Gctf","url":"https://www.mrc-lmb.cam.ac.uk/kzhang/Gctf/","desc":"A GPU-accelerated CTF estimation tool with high speed and precision"},{"name":"Warp","url":"http://www.warpem.com/","desc":"A fast preprocessing tool integrating motion correction, CTF estimation, and particle picking"}]},{"title":"Particle Picking","items":[{"name":"crYOLO","url":"https://cryolo.readthedocs.io/","desc":"A deep learning-based particle picking tool supporting real-time picking"},{"name":"Topaz","url":"https://topaz.readthedocs.io/","desc":"A CNN-based particle picking tool with excellent performance on low signal-to-noise ratio samples"}]},{"title":"Modeling and Visualization","items":[{"name":"PyMOL","url":"https://pymol.org/","desc":"A molecular visualization software, the standard tool for scripted rendering and structure analysis"},{"name":"UCSF ChimeraX","url":"https://www.cgl.ucsf.edu/chimerax/","desc":"A modern molecular visualization software with outstanding density map processing capabilities"},{"name":"Coot","url":"https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/","desc":"A manual model building tool for interactive adjustment of density maps and models"},{"name":"phenix","url":"https://phenix-online.org/","desc":"A structure refinement and validation suite supporting cryo-EM and crystallography"},{"name":"ISOLDE","url":"https://isolde.cimr.cam.ac.uk/","desc":"An interactive real-time molecular dynamics modeling tool based on ChimeraX"},{"name":"ModelAngelo","url":"https://github.com/3dem/model-angelo","desc":"A deep learning-based automatic atomic modeling tool for density maps"}]}]},{"title":"Database Resources","children":[{"title":"Structures and Density Maps","items":[{"name":"EMDB","url":"https://www.ebi.ac.uk/emdb/","desc":"The Electron Microscopy Data Bank, storing and distributing 3D density maps from cryo-EM"},{"name":"EMPIAR","url":"https://www.ebi.ac.uk/empiar/","desc":"The Electron Microscopy Public Image Archive, providing raw movies and particle data"},{"name":"RCSB PDB","url":"https://www.rcsb.org/","desc":"The Protein Data Bank, containing biological macromolecular structures determined by X-ray, NMR, and cryo-EM"},{"name":"wwPDB","url":"https://www.wwpdb.org/","desc":"The worldwide PDB coordination body, uniformly managing structural data standards"},{"name":"PDBe","url":"https://www.ebi.ac.uk/pdbe/","desc":"The European PDB node, providing structural data and visualization tools"},{"name":"AlphaFold DB","url":"https://alphafold.ebi.ac.uk/","desc":"The AlphaFold protein structure database, covering hundreds of millions of proteins"}]}]}]}]`),ng=Object.freeze(Object.defineProperty({__proto__:null,default:sg},Symbol.toStringTag,{value:"Module"})),ag=JSON.parse('[{"title":"数据分析","children":[{"title":"数据探索","children":[{"title":"R 语言","items":[{"name":"tidyverse","url":"https://www.tidyverse.org/","desc":"R 语言中用于数据科学的一整套包，涵盖数据导入、清洗、转换和可视化"},{"name":"data.table","url":"https://rdatatable.gitlab.io/data.table/","desc":"高性能数据操作包，百万行级数据框处理速度极快"},{"name":"readr","url":"https://readr.tidyverse.org/","desc":"快速读取 CSV/TSV 等表格文件的工具，自动推断列类型"}]},{"title":"Python","items":[{"name":"NumPy","url":"https://numpy.org/","desc":"Python 科学计算基础库，提供多维数组与向量化运算"},{"name":"Pandas","url":"https://pandas.pydata.org/","desc":"基于 Python 的强大数据分析库，提供高效的数据操作和处理功能"},{"name":"Polars","url":"https://www.pola.rs/","desc":"基于 Apache Arrow 的高性能数据处理库，支持多语言接口"},{"name":"Jupyter","url":"https://jupyter.org/","desc":"交互式笔记本环境，支持代码、文本与可视化混排"}]}]},{"title":"数据可视化","children":[{"title":"R 语言","items":[{"name":"ggplot2","url":"https://ggplot2.tidyverse.org/","desc":"R 语言中功能强大的数据可视化包，基于语法图形学"},{"name":"plotly（R）","url":"https://plotly.com/r/","desc":"交互式图表库的 R 接口，可生成网页端可交互图形"},{"name":"patchwork","url":"https://patchwork.data-imaginist.com/","desc":"多图拼接工具，用 + 和 / 符号轻松组合 ggplot 图形"},{"name":"RColorBrewer","url":"https://cran.r-project.org/web/packages/RColorBrewer/","desc":"出版级配色方案，提供精心设计的离散色板"}]},{"title":"Python","items":[{"name":"Matplotlib","url":"https://matplotlib.org/","desc":"Python 中广泛使用的绘图库，支持多种图表类型"},{"name":"Seaborn","url":"https://seaborn.pydata.org/","desc":"基于 Matplotlib 的统计可视化库，内置美观主题与统计图"},{"name":"Plotly","url":"https://plotly.com/python/","desc":"交互式可视化库，支持散点、3D、地图等多种图表"},{"name":"Altair","url":"https://altair-viz.github.io/","desc":"基于 Vega-Lite 的声明式统计可视化库"}]}]},{"title":"统计分析","children":[{"title":"R 语言","items":[{"name":"R Stats","url":"https://www.r-project.org/","desc":"R 语言内置的统计分析功能，涵盖广泛的统计方法和模型"},{"name":"rstatix","url":"https://rpkgs.datanovia.com/rstatix/","desc":"管道友好的统计检验封装（t 检验、方差分析、秩和检验等）"},{"name":"broom","url":"https://broom.tidymodels.org/","desc":"把统计模型输出整理为整洁数据框的工具"},{"name":"lme4","url":"https://cran.r-project.org/web/packages/lme4/","desc":"线性与非线性混合效应模型，适用于重复测量数据"}]},{"title":"Python","items":[{"name":"SciPy Stats","url":"https://scipy.org/","desc":"Python 科学计算库 SciPy 中的统计模块，提供多种统计分布和检验方法"},{"name":"statsmodels","url":"https://www.statsmodels.org/","desc":"统计建模库，支持回归、时间序列、假设检验等"},{"name":"pingouin","url":"https://pingouin-stats.org/","desc":"简洁易用的统计检验库，覆盖常用参数与非参数检验"}]}]},{"title":"机器学习","children":[{"title":"Python","items":[{"name":"scikit-learn","url":"https://scikit-learn.org/","desc":"基于 Python 的机器学习库，提供丰富的算法和工具"},{"name":"PyTorch","url":"https://pytorch.org/","desc":"深度学习框架，支持动态计算图和高效的张量运算"},{"name":"TensorFlow","url":"https://www.tensorflow.org/","desc":"深度学习框架，生态成熟，支持生产级部署"},{"name":"XGBoost","url":"https://xgboost.readthedocs.io/","desc":"梯度提升树库，表格数据竞赛与工程的首选"},{"name":"LightGBM","url":"https://lightgbm.readthedocs.io/","desc":"微软出品的梯度提升框架，训练速度快、内存占用低"}]},{"title":"R 语言","items":[{"name":"tidymodels","url":"https://www.tidymodels.org/","desc":"R 语言统一的建模框架，覆盖数据预处理、建模与评估"}]}]}]},{"title":"组学","children":[{"title":"数据平台","children":[{"title":"序列数据库","items":[{"name":"NCBI","url":"https://www.ncbi.nlm.nih.gov/","desc":"美国国家生物技术信息中心，提供 GenBank 等重要数据库"},{"name":"EMBL-EBI","url":"https://www.ebi.ac.uk/","desc":"欧洲分子生物学实验室，提供多种生物信息学资源和工具"},{"name":"DDBJ","url":"https://www.ddbj.nig.ac.jp/","desc":"日本 DNA 数据库，提供核酸、蛋白序列数据的存储和访问"},{"name":"UniProt","url":"https://www.uniprot.org/","desc":"最全面的蛋白质序列与功能注释数据库"},{"name":"Ensembl","url":"https://www.ensembl.org/","desc":"脊椎动物基因组注释与比较基因组学数据库"},{"name":"UCSC Genome Browser","url":"https://genome.ucsc.edu/","desc":"基因组浏览器，支持多物种基因组可视化与自定义轨道"}]},{"title":"测序数据","items":[{"name":"SRA","url":"https://www.ncbi.nlm.nih.gov/sra","desc":"NCBI 原始测序数据仓库，存储高通量测序原始数据"},{"name":"GEO","url":"https://www.ncbi.nlm.nih.gov/geo/","desc":"基因表达数据库，收录芯片与测序表达数据"},{"name":"ENA","url":"https://www.ebi.ac.uk/ena","desc":"欧洲核酸档案库，欧洲的原始测序数据存储中心"}]},{"title":"蛋白互作与通路","items":[{"name":"STRING","url":"https://string-db.org/","desc":"蛋白-蛋白相互作用数据库，支持多种物种"},{"name":"KEGG","url":"https://www.kegg.jp/","desc":"通路、基因与化合物综合数据库"},{"name":"GO","url":"http://geneontology.org/","desc":"基因本体论，提供标准化的基因功能注释体系"}]}]},{"title":"分析软件","children":[{"title":"命令行","items":[{"name":"FastQC","url":"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/","desc":"测序数据质量评估工具，生成可视化质控报告"},{"name":"Trimmomatic","url":"http://www.usadellab.org/cms/?page=trimmomatic","desc":"测序读段修剪工具，去除接头与低质量碱基"},{"name":"cutadapt","url":"https://cutadapt.readthedocs.io/","desc":"命令行接头去除工具，支持多种测序平台"},{"name":"BWA","url":"https://github.com/lh3/bwa","desc":"短读段比对工具，基因组比对的标准选择"},{"name":"STAR","url":"https://github.com/alexdobin/STAR","desc":"RNA-seq 剪接感知比对器，速度极快"},{"name":"HISAT2","url":"https://daehwankimlab.github.io/hisat2/","desc":"基于图索引的 RNA-seq 比对工具"},{"name":"SAMtools","url":"https://www.htslib.org/","desc":"SAM/BAM 文件处理工具集，比对后处理的标准"},{"name":"bcftools","url":"https://samtools.github.io/bcftools/","desc":"VCF/BCF 变异文件处理工具，支持过滤与统计"},{"name":"GATK","url":"https://gatk.broadinstitute.org/","desc":"基因组分析工具包，变异检测的行业标准"},{"name":"MultiQC","url":"https://multiqc.info/","desc":"把多个质控报告汇总为一个交互式 HTML 报告"}]},{"title":"Python","items":[{"name":"Biopython","url":"https://biopython.org/","desc":"生物信息学核心 Python 库，涵盖序列、结构与文件解析"},{"name":"scikit-bio","url":"https://scikit-bio.org/","desc":"用于生物信息学的 Python 库，提供多种生物数据处理和分析工具"},{"name":"scanpy","url":"https://scanpy.readthedocs.io/","desc":"单细胞转录组分析 Python 库，与 AnnData 深度集成"}]},{"title":"R 语言","items":[{"name":"Bioconductor","url":"https://bioconductor.org/","desc":"基于 R 语言的生物信息学软件包集合，涵盖基因组学、转录组学等领域"},{"name":"Seurat","url":"https://satijalab.org/seurat/","desc":"单细胞转录组分析主流 R 包，提供完整分析流程"},{"name":"DESeq2","url":"https://bioconductor.org/packages/release/bioc/html/DESeq2.html","desc":"基于负二项模型的差异表达分析标准工具"},{"name":"edgeR","url":"https://bioconductor.org/packages/release/bioc/html/edgeR.html","desc":"基于经验贝叶斯的差异表达分析工具"},{"name":"limma","url":"https://bioconductor.org/packages/release/bioc/html/limma.html","desc":"线性模型差异分析工具，芯片与测序数据通用"}]}]}]},{"title":"电镜结构解析","children":[{"title":"软件工具","children":[{"title":"三维重构","items":[{"name":"RELION","url":"https://www3.mrc-lmb.cam.ac.uk/groups/scheres/index.php","desc":"基于贝叶斯方法的冷冻电镜三维重构软件，广泛应用于高分辨率结构解析"},{"name":"cryoSPARC","url":"https://cryoem.slac.stanford.edu/","desc":"高性能的冷冻电镜数据处理平台，提供直观的用户界面"},{"name":"cisTEM","url":"https://cistem.org/","desc":"免费开源的冷冻电镜处理软件，支持完整的单颗粒流程"},{"name":"SPHIRE","url":"https://sphire.mpg.de/","desc":"模块化冷冻电镜处理套件，支持高分辨率重构"}]},{"title":"预处理","items":[{"name":"MotionCor2","url":"https://emcore.ucsf.edu/ucsf-software","desc":"逐帧运动校正工具，消除电子束导致的样品漂移"},{"name":"CTFFIND4","url":"https://grigoriefflab.umassmed.edu/ctffind4","desc":"CTF 估计经典工具，评估欠焦与像散"},{"name":"Gctf","url":"https://www.mrc-lmb.cam.ac.uk/kzhang/Gctf/","desc":"GPU 加速的 CTF 估计工具，速度快精度高"},{"name":"Warp","url":"http://www.warpem.com/","desc":"集运动校正、CTF 与颗粒挑选于一体的快速预处理工具"}]},{"title":"颗粒挑选","items":[{"name":"crYOLO","url":"https://cryolo.readthedocs.io/","desc":"基于深度学习的颗粒挑选工具，支持实时挑选"},{"name":"Topaz","url":"https://topaz.readthedocs.io/","desc":"基于 CNN 的颗粒挑选工具，对低信噪比样品表现优异"}]},{"title":"建模与可视化","items":[{"name":"PyMOL","url":"https://pymol.org/","desc":"分子可视化软件，脚本化渲染与结构分析的标准工具"},{"name":"UCSF ChimeraX","url":"https://www.cgl.ucsf.edu/chimerax/","desc":"现代分子可视化软件，密度图处理能力突出"},{"name":"Coot","url":"https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/","desc":"手动模型修正工具，密度图与模型交互调整"},{"name":"phenix","url":"https://phenix-online.org/","desc":"结构精修与验证套件，支持冷冻电镜与晶体学"},{"name":"ISOLDE","url":"https://isolde.cimr.cam.ac.uk/","desc":"交互式实时分子动力学建模工具，基于 ChimeraX"},{"name":"ModelAngelo","url":"https://github.com/3dem/model-angelo","desc":"基于深度学习的密度图自动原子建模工具"}]}]},{"title":"数据库资源","children":[{"title":"结构与密度图","items":[{"name":"EMDB","url":"https://www.ebi.ac.uk/emdb/","desc":"电子显微镜数据银行，存储和分发冷冻电镜三维密度图"},{"name":"EMPIAR","url":"https://www.ebi.ac.uk/empiar/","desc":"冷冻电镜原始数据档案，提供原始电影与颗粒数据"},{"name":"RCSB PDB","url":"https://www.rcsb.org/","desc":"蛋白质数据银行，包含通过X射线、NMR和冷冻电镜确定的生物大分子结构"},{"name":"wwPDB","url":"https://www.wwpdb.org/","desc":"全球 PDB 协调机构，统一管理结构数据标准"},{"name":"PDBe","url":"https://www.ebi.ac.uk/pdbe/","desc":"欧洲 PDB 节点，提供结构数据与可视化工具"},{"name":"AlphaFold DB","url":"https://alphafold.ebi.ac.uk/","desc":"AlphaFold 预测蛋白质结构数据库，覆盖数亿蛋白质"}]}]}]}]'),eg=Object.freeze(Object.defineProperty({__proto__:null,default:ag},Symbol.toStringTag,{value:"Module"})),tg=[{name:"Advanced",count:1},{name:"AlphaFold",count:1},{name:"Atomic Modeling",count:1},{name:"AutoDock Vina",count:1},{name:"Bash",count:1},{name:"ChimeraX",count:1},{name:"Command Line",count:1},{name:"Computational Biology",count:2},{name:"Coot",count:1},{name:"Cryo-Electron Microscopy",count:2},{name:"cryo-EM",count:2},{name:"Data Processing",count:2},{name:"Data Visualization",count:1},{name:"dplyr",count:1},{name:"Getting Started",count:2},{name:"ggplot2",count:1},{name:"Linux",count:1},{name:"Molecular Docking",count:1},{name:"NumPy",count:1},{name:"Pandas",count:1},{name:"phenix",count:1},{name:"Protein Design",count:1},{name:"PyMOL",count:1},{name:"Python",count:3},{name:"R Language",count:3},{name:"RELION",count:1},{name:"Review",count:1},{name:"RFdiffusion",count:1},{name:"Rosetta",count:1},{name:"Scripting",count:1},{name:"Shell",count:1},{name:"Structure Refinement",count:1},{name:"Structure Visualization",count:1},{name:"tidyverse",count:1},{name:"Tutorial",count:10},{name:"Virtual Screening",count:1}],lg=Object.freeze(Object.defineProperty({__proto__:null,default:tg},Symbol.toStringTag,{value:"Module"})),og=[{name:"蛋白质设计",count:1},{name:"分子对接",count:1},{name:"计算生物学",count:2},{name:"脚本编程",count:1},{name:"教程",count:10},{name:"结构精修",count:1},{name:"结构可视化",count:1},{name:"进阶",count:1},{name:"冷冻电镜",count:2},{name:"命令行",count:1},{name:"入门",count:2},{name:"数据处理",count:2},{name:"数据可视化",count:1},{name:"虚拟筛选",count:1},{name:"原子建模",count:1},{name:"综述",count:1},{name:"AlphaFold",count:1},{name:"AutoDock Vina",count:1},{name:"Bash",count:1},{name:"ChimeraX",count:1},{name:"Coot",count:1},{name:"cryo-EM",count:2},{name:"dplyr",count:1},{name:"ggplot2",count:1},{name:"Linux",count:1},{name:"NumPy",count:1},{name:"Pandas",count:1},{name:"phenix",count:1},{name:"PyMOL",count:1},{name:"Python",count:3},{name:"R语言",count:3},{name:"RELION",count:1},{name:"RFdiffusion",count:1},{name:"Rosetta",count:1},{name:"Shell",count:1},{name:"tidyverse",count:1}],pg=Object.freeze(Object.defineProperty({__proto__:null,default:og},Symbol.toStringTag,{value:"Module"})),cg=[{name:"StructureOfProteinDemo",title:"Protein Structure Determination Demo Project",desc:"Protein structure placeholder demo",date:"2026-04-12",tags:["Protein","Structure","Cryo-EM"],doi:"10.1234/demo-structure.2026",url:"https://www.demo-structure.org",status:"completed",journal:"Journal of Demo Science",year:2026,authors:["Demo Author A","Demo Author B"]}],rg=Object.freeze(Object.defineProperty({__proto__:null,default:cg},Symbol.toStringTag,{value:"Module"})),ig=[{name:"StructureOfProteinDemo",title:"蛋白质结构解析演示课题",desc:"蛋白质结构占位demo",date:"2026-04-12",tags:["Protein","Structure","Cryo-EM"],doi:"10.1234/demo-structure.2026",url:"https://www.demo-structure.org",status:"completed",journal:"Journal of Demo Science",year:2026,authors:["Demo Author A","Demo Author B"]}],ug=Object.freeze(Object.defineProperty({__proto__:null,default:ig},Symbol.toStringTag,{value:"Module"})),hg=`<h1>Mainstream Tools for Protein Design</h1>
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
<p>The next article will introduce mainstream tools for virtual screening: a complete toolchain for docking small molecules to targets.</p>`,dg=`<h1>蛋白质设计主流工具</h1>
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
<p>下一篇将介绍虚拟筛选的主流工具：把小分子与靶点对接的完整工具链。</p>`,mg=`<h1>Mainstream Tools for Virtual Screening</h1>
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
<p>This completes the computational biology direction: protein design tools + virtual screening tools, with two main threads forming a closed loop of "design-screen" computational drug discovery.</p>`,jg=`<h1>虚拟筛选主流工具</h1>
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
<p>至此计算生物学方向完成：蛋白质设计工具 + 虚拟筛选工具，两条主线构成"设计-筛选"的计算药物发现闭环。</p>`,gg=`<h1>Bash Programming: Variables, Conditionals, Loops, and Practical Scripts</h1>
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
<p>At this point, the programming track tutorials are complete: Python ×3, R ×3, Linux ×1, Bash ×1, from zero to hands-on practice.</p>`,fg=`<h1>Bash 编程：变量、条件、循环与实用脚本</h1>
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
<p>至此编程方向教程完成：Python ×3、R ×3、Linux ×1、Bash ×1，从零到实战。</p>`,bg=`<h1>Linux Command Line Basics: File System, Permissions, and Text Processing</h1>
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
<p>The next article will introduce Bash programming: organizing commands into reusable scripts.</p>`,_g=`<h1>Linux 命令行基础：文件系统、权限与文本处理</h1>
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
<p>下一篇将介绍 Bash 编程：把命令组织成可复用的脚本。</p>`,yg=`<h1>Python Advanced: Functions, Classes, and Modules</h1>
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
<p>The next article will cover Python data processing in practice: file I/O, regular expressions, and NumPy/Pandas.</p>`,vg=`<h1>Python 进阶：函数、类与模块</h1>
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
<p>下一篇将介绍 Python 数据处理实战：文件 IO、正则表达式与 NumPy/Pandas。</p>`,wg=`<h1>Introduction to Python Programming: Environment, Syntax, and Data Types</h1>
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
<p>The next article will introduce functions, classes, and modules, moving into real engineering-style programming.</p>`,xg=`<h1>Python 编程入门：环境、语法与数据类型</h1>
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
<p>下一篇将介绍函数、类与模块，进入真正的工程化编程。</p>`,Cg=`<h1>Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas</h1>
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
<p>At this point, the Python tutorial trilogy is complete: Beginner → Intermediate → Data Processing.</p>`,kg=`<h1>Python 数据处理实战：文件 IO、正则与 NumPy/Pandas</h1>
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
<p>至此 Python 教程三部曲完成：入门 → 进阶 → 数据处理。</p>`,Sg=`<h1>R Language Primer: Data Structures, Vectorization, and Functions</h1>
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
<p>The next article will introduce tidyverse: using <code>dplyr</code> for elegant data manipulation.</p>`,Pg=`<h1>R 语言入门：数据结构、向量化与函数</h1>
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
<p>下一篇将介绍 tidyverse：用 <code>dplyr</code> 优雅地完成数据操作。</p>`,Ag=`<h1>R ggplot2 Data Visualization: Grammar of Graphics and Common Charts</h1>
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
<p>With this, the R tutorial trilogy is complete: Introduction → tidyverse → ggplot2.</p>`,Tg=`<h1>R ggplot2 数据可视化：图层语法与常用图表</h1>
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
<p>至此 R 教程三部曲完成：入门 → tidyverse → ggplot2。</p>`,Rg=`<h1>R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning</h1>
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
<p>The next article will introduce ggplot2: creating publication-quality graphics using the grammar of layers.</p>`,Eg=`<h1>R tidyverse 数据操作：dplyr 管道与数据清洗</h1>
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
<p>下一篇将介绍 ggplot2：用图层语法绘制出版级图表。</p>`,Mg=`<h1>Review of Cryo-EM Technology</h1>
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
<p>The past decade of cryo-EM is a model of synergy among engineering, physics, and computational biology. In the next decade, it will join forces with AI to redefine how we understand life.</p>`,Dg=`<h1>冷冻电镜技术综述</h1>
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
<p>冷冻电镜的十年，是工程、物理与计算生物学协同的典范。下一个十年，它将与 AI 一起重新定义我们理解生命的方式。</p>`,Lg=`<h1>Cryo-EM Single Particle Analysis: Full Data Processing Workflow</h1>
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
<p>Future articles will cover structure visualization (PyMOL/ChimeraX) and atomic modeling (Coot/phenix), turning density maps into atomic models.</p>`,Og=`<h1>冷冻电镜单颗粒分析：数据处理全流程</h1>
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
<p>后续文章将介绍结构可视化（PyMOL/ChimeraX）与原子建模（Coot/phenix），把密度图变成原子模型。</p>`,Fg=`<h1>Atomic Modeling and Refinement</h1>
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
<p>This completes the structural biology pipeline: data processing → technical review → visualization → atomic modeling, forming a complete loop from experiment to model.</p>`,$g=`<h1>原子建模与精修</h1>
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
<p>至此结构生物学方向完成：数据处理流程 → 技术综述 → 可视化 → 原子建模，构成从实验到模型的完整闭环。</p>`,Ig=`<h1>Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice</h1>
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
<p>The next article will cover atomic modeling: how to build and refine atomic models from density maps.</p>`,Ng=`<h1>生物大分子结构可视化：PyMOL 与 ChimeraX 实战</h1>
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
<p>下一篇将介绍原子建模：如何从密度图搭建并精修原子模型。</p>`,Il=()=>{const s=globalThis.__GBLOG_LOCALE__;return s||(typeof window<"u"?localStorage.getItem("locale"):null)||"zh-CN"},Bg=(s,n=".json")=>{const e=Il()==="en-US";return s.endsWith("-en")?`${s}${n}`:e?`${s}-en${n}`:`${s}${n}`},Se=Object.assign({"../content/about-en.json":Oj,"../content/about.json":$j,"../content/categories-en.json":Nj,"../content/categories.json":zj,"../content/notes-en.json":Hj,"../content/notes.json":Gj,"../content/posts-en.json":Wj,"../content/posts.json":Yj,"../content/projects-en.json":Jj,"../content/projects.json":Zj,"../content/resources-en.json":ng,"../content/resources.json":eg,"../content/tags-en.json":lg,"../content/tags.json":pg,"../content/topics-en.json":rg,"../content/topics.json":ug}),Pp=Object.assign({"../content/html/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en.html":hg,"../content/html/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools.html":dg,"../content/html/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools-en.html":mg,"../content/html/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools.html":jg,"../content/html/notes/Programming/bash/bash-scripting/bash-scripting-en.html":gg,"../content/html/notes/Programming/bash/bash-scripting/bash-scripting.html":fg,"../content/html/notes/Programming/linux/linux-basics/linux-basics-en.html":bg,"../content/html/notes/Programming/linux/linux-basics/linux-basics.html":_g,"../content/html/notes/Programming/python/python-advanced/python-advanced-en.html":yg,"../content/html/notes/Programming/python/python-advanced/python-advanced.html":vg,"../content/html/notes/Programming/python/python-basics/python-basics-en.html":wg,"../content/html/notes/Programming/python/python-basics/python-basics.html":xg,"../content/html/notes/Programming/python/python-data/python-data-en.html":Cg,"../content/html/notes/Programming/python/python-data/python-data.html":kg,"../content/html/notes/Programming/r/r-basics/r-basics-en.html":Sg,"../content/html/notes/Programming/r/r-basics/r-basics.html":Pg,"../content/html/notes/Programming/r/r-ggplot2/r-ggplot2-en.html":Ag,"../content/html/notes/Programming/r/r-ggplot2/r-ggplot2.html":Tg,"../content/html/notes/Programming/r/r-tidyverse/r-tidyverse-en.html":Rg,"../content/html/notes/Programming/r/r-tidyverse/r-tidyverse.html":Eg,"../content/html/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en.html":Mg,"../content/html/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview.html":Dg,"../content/html/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow-en.html":Lg,"../content/html/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow.html":Og,"../content/html/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling-en.html":Fg,"../content/html/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling.html":$g,"../content/html/notes/StructuralBiology/visualization/structure-visualization/structure-visualization-en.html":Ig,"../content/html/notes/StructuralBiology/visualization/structure-visualization/structure-visualization.html":Ng}),Da=s=>{const n=Bg(s,".json"),a=Object.keys(Se).find(e=>e.includes(n));if(a)return Se[a].default||{};if(console.error(`Failed to load JSON content: ${n}`),Il()==="en-US"){const e=Object.keys(Se).find(t=>t.includes(`${s}.json`));if(e)return Se[e].default||{}}return{}},zg=s=>{const a=Il()==="en-US",e=s.replace(/\.md$/i,""),t=[];e.endsWith("-en")?(t.push(e),a&&t.push(e.replace(/-en$/,""))):(t.push(a?`${e}-en`:e),t.push(e));for(const l of t){const o=Object.keys(Pp).find(p=>p.includes(`/content/html/${l}.html`));if(o)return Pp[o]}return`<h1>Content Not Available</h1>
<p>The requested content could not be loaded.</p>`},ut=()=>Da("categories"),qg=()=>Da("posts"),Hg=()=>Da("notes"),Vg=()=>Da("tags"),Gg=()=>Da("about"),Ug=()=>Da("resources"),Wg={class:"navigation-tree"},Kg={class:"tree-label"},Yg={class:"article-list article-list-root"},Xg={key:0,class:"article-item"},Jg={class:"directory-node"},Qg={class:"tree-item tree-item--folder"},Zg={key:0,class:"article-list sub-list files-level"},sf={key:1,class:"article-list sub-list"},nf={class:"directory-node"},af={class:"tree-item tree-item--folder"},ef={key:0,class:"article-list sub-list files-level"},tf=zs({__name:"NavigationTree",setup(s){const{locale:n}=Xs(),a=Ma(),e=ns([]),t=ns(""),l=ns([]);function o(){try{l.value=ut()||[]}catch(i){console.error("Failed to load category data:",i),l.value=[]}}function p(i){return $n(`/article/${i.replace(/\.md$/,"")}`)}function c(i){const d=t.value.replace(/\.md$/,""),h=i.replace(/\.md$/,""),j=d.replace(/-en$/,""),m=h.replace(/-en$/,"");return j===m}function u(){const i=t.value;if(!i){e.value=[];return}const d=i.split("/").filter(Boolean);if(d.length<2){e.value=[];return}const h=d[0],j=d[1];let m=null;if(Array.isArray(l.value)){s:for(const x of l.value)if(Array.isArray(x.items))for(const f of x.items){if(f?.name!==j)continue;let k=!1;if(Array.isArray(f.articles)&&(k=f.articles.some(w=>typeof w?.articleUrl=="string"&&w.articleUrl.includes(`/article/${h}/`))),!k&&Array.isArray(f.categories)&&(k=f.categories.some(w=>Array.isArray(w.articles)&&w.articles.some(y=>typeof y?.articleUrl=="string"&&y.articleUrl.includes(`/article/${h}/${j}/`)))),k){m=f;break s}}}if(!m){e.value=[];return}const P=[],C=[],b=(x,f)=>{const k=String(f).replace(/^\/+/,"").split("/"),w=k[0]==="article"?1:0,y=k[w],M=k[w+1],S=k.slice(w+2),F=`${y}/${M}/${S.join("/")}`;return{title:x,path:`${F}.md`}};Array.isArray(m.articles)&&m.articles.forEach(x=>{x?.articleUrl&&P.push(b(x.title,x.articleUrl))}),Array.isArray(m.categories)&&m.categories.forEach(x=>{const f=[];Array.isArray(x.articles)&&x.articles.forEach(k=>{k?.articleUrl&&f.push(b(k.title,k.articleUrl))}),f.length&&C.push({name:x.title||x.key,type:"directory",files:f})}),e.value=[{name:m.title||m.name||j,files:P,children:C}]}function r(){const i=a.params.path;t.value=Array.isArray(i)?i.join("/"):i||"",u()}return Bs(()=>a.params.path,r,{immediate:!0}),Bs(n,()=>{o(),r()}),dn(()=>{o(),u()}),(i,d)=>{const h=de("router-link");return E(),L("div",Wg,[(E(!0),L(js,null,Ts(e.value,j=>(E(),L("div",{key:j.name,class:"tree-group"},[v("div",Kg,G(j.name),1),v("ul",Yg,[(E(!0),L(js,null,Ts(j.children,m=>(E(),L(js,{key:m.name},[m.files&&m.files.length?(E(),L("li",Xg,[v("div",Jg,[v("div",Qg,G(m.name),1),m.files&&m.files.length?(E(),L("ul",Zg,[(E(!0),L(js,null,Ts(m.files,P=>(E(),L("li",{key:P.path,class:"article-item"},[bs(h,{to:p(P.path),class:Ns(["tree-item tree-item--child",{"tree-item--active":c(P.path)}])},{default:hn(()=>[Ss(G(P.title),1)]),_:2},1032,["to","class"])]))),128))])):rs("",!0),m.children&&m.children.length?(E(),L("ul",sf,[(E(!0),L(js,null,Ts(m.children,P=>(E(),L("li",{key:P.name,class:"article-item"},[v("div",nf,[v("div",af,G(P.name),1),P.files&&P.files.length?(E(),L("ul",ef,[(E(!0),L(js,null,Ts(P.files,C=>(E(),L("li",{key:C.path,class:"article-item"},[bs(h,{to:p(C.path),class:Ns(["tree-item tree-item--child",{"tree-item--active":c(C.path)}])},{default:hn(()=>[Ss(G(C.title),1)]),_:2},1032,["to","class"])]))),128))])):rs("",!0)])]))),128))])):rs("",!0)])])):rs("",!0)],64))),128)),(E(!0),L(js,null,Ts(j.files,m=>(E(),L("li",{key:m.path,class:"article-item"},[bs(h,{to:p(m.path),class:Ns(["tree-item tree-item--child",{"tree-item--active":c(m.path)}])},{default:hn(()=>[Ss(G(m.title),1)]),_:2},1032,["to","class"])]))),128))])]))),128))])}}}),mn=(s,n)=>{const a=s.__vccOpts||s;for(const[e,t]of n)a[e]=t;return a},qr=mn(tf,[["__scopeId","data-v-7ae1de2a"]]),lf={class:"app-header"},of={class:"container px-0"},pf={class:"navbar navbar-expand-lg app-nav"},cf={class:"container-fluid d-flex app-nav__inner"},rf={class:"app-nav__wordmark"},uf={class:"d-flex d-lg-none ms-auto app-nav__actions app-nav__actions--mobile"},hf=["aria-label"],df=["aria-label"],mf={key:0,class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},jf={key:1,class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},gf=["aria-label"],ff=["aria-label"],bf={class:"navbar-nav mb-2 mb-lg-0 app-nav__links"},_f={key:0,class:"app-nav__link-divider","aria-hidden":"true"},yf={class:"d-none d-lg-flex ms-auto app-nav__actions"},vf=["aria-label"],wf=["aria-label"],xf={key:0,class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},Cf={key:1,class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},kf=["aria-label"],Sf={class:"offcanvas-panel"},Pf={class:"offcanvas-section"},Af={class:"offcanvas-head"},Tf={class:"offcanvas-card"},Rf={class:"list-unstyled m-0"},Ef={key:0,class:"offcanvas-section"},Mf={class:"offcanvas-head"},Df=zs({__name:"AppHeader",emits:["open-search"],setup(s,{emit:n}){const a=Ma(),e=fe(),{t}=Xs(),l=$l(),o=ss(()=>l.theme==="dark"?!0:l.theme==="light"?!1:typeof window<"u"&&window.matchMedia("(prefers-color-scheme: dark)").matches),p=ns(!1),c=ns(!1),u=n,r=w=>{w.currentTarget?.blur()},i=ss(()=>[{text:t("categories"),href:$n("/category")},{text:t("resources"),href:$n("/resource")},{text:t("about"),href:$n("/about")}]),d=w=>w===a.path?!0:w!==$n("/")&&a.path.startsWith(w),h=ss(()=>a.path.includes("/article/")),j=()=>{l.toggleTheme()},m=()=>{const w=a.path.match(/^\/(zh|en)/),y=w&&w[1]==="en"?"/zh":"/en",M=w?a.path.slice(w[0].length):a.path;e.push({path:`${y}${M}`,query:a.query})},P=()=>{window.innerWidth<992?x():p.value=!p.value},C=()=>{const w=document.documentElement,y=document.body;w&&(w.style.overflow="hidden",w.style.overscrollBehavior="contain"),y&&(y.style.overflow="hidden",y.style.overscrollBehavior="contain")},b=()=>{const w=document.documentElement,y=document.body;w&&(w.style.overflow="",w.style.overscrollBehavior=""),y&&(y.style.overflow="",y.style.overscrollBehavior="")},x=()=>{C(),c.value=!0},f=()=>{c.value=!1,b()},k=w=>{const y=w.target;y&&y.closest("a")&&f()};return dn(()=>{l.initTheme(),l.initLocale()}),Tn(()=>{b()}),(w,y)=>{const M=de("RouterLink");return E(),L(js,null,[v("header",lf,[v("div",of,[v("nav",pf,[v("div",cf,[bs(M,{class:"navbar-brand app-nav__brand",to:ls($n)("/"),onClick:y[0]||(y[0]=S=>p.value=!1)},{default:hn(()=>[v("span",rf,[Ss(G(ls(en).author),1),y[4]||(y[4]=v("span",{class:"app-nav__apos"},"’",-1)),y[5]||(y[5]=Ss("s blog",-1))])]),_:1},8,["to"]),v("div",uf,[v("button",{class:"icon-btn",onClick:y[1]||(y[1]=S=>u("open-search")),onFocus:r,"aria-label":ls(t)("search")},[...y[6]||(y[6]=[v("svg",{class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[v("circle",{cx:"11",cy:"11",r:"8"}),v("line",{x1:"21",y1:"21",x2:"16.65",y2:"16.65"})],-1)])],40,hf),v("button",{class:"icon-btn",onClick:j,onFocus:r,"aria-label":ls(t)("theme")},[o.value?(E(),L("svg",mf,[...y[7]||(y[7]=[Ce('<circle cx="12" cy="12" r="5" data-v-b2e68a7c></circle><line x1="12" y1="1" x2="12" y2="3" data-v-b2e68a7c></line><line x1="12" y1="21" x2="12" y2="23" data-v-b2e68a7c></line><line x1="4.22" y1="4.22" x2="5.64" y2="5.64" data-v-b2e68a7c></line><line x1="18.36" y1="18.36" x2="19.78" y2="19.78" data-v-b2e68a7c></line><line x1="1" y1="12" x2="3" y2="12" data-v-b2e68a7c></line><line x1="21" y1="12" x2="23" y2="12" data-v-b2e68a7c></line><line x1="4.22" y1="19.78" x2="5.64" y2="18.36" data-v-b2e68a7c></line><line x1="18.36" y1="5.64" x2="19.78" y2="4.22" data-v-b2e68a7c></line>',9)])])):(E(),L("svg",jf,[...y[8]||(y[8]=[v("path",{d:"M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"},null,-1)])]))],40,df),v("button",{class:"icon-btn",onClick:m,onFocus:r,"aria-label":ls(t)("language")},[...y[9]||(y[9]=[Ce('<svg class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.75" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" data-v-b2e68a7c><path d="m5 8 6 6" data-v-b2e68a7c></path><path d="m4 14 6-6 2-3" data-v-b2e68a7c></path><path d="M2 5h12" data-v-b2e68a7c></path><path d="M7 2h1" data-v-b2e68a7c></path><path d="m22 22-5-10-5 10" data-v-b2e68a7c></path><path d="M14 18h6" data-v-b2e68a7c></path></svg>',1)])],40,gf),v("button",{class:"icon-btn",onClick:P,onFocus:r,"aria-label":ls(t)("menu")},[...y[10]||(y[10]=[v("svg",{class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[v("line",{x1:"3",y1:"12",x2:"21",y2:"12"}),v("line",{x1:"3",y1:"6",x2:"21",y2:"6"}),v("line",{x1:"3",y1:"18",x2:"21",y2:"18"})],-1)])],40,ff)]),v("div",{class:Ns(["navbar-collapse collapse",{show:p.value}])},[v("ul",bf,[i.value.length?(E(),L("li",_f)):rs("",!0),(E(!0),L(js,null,Ts(i.value,S=>(E(),L("li",{class:"nav-item",key:S.text},[bs(M,{class:Ns(["nav-link app-nav__link",{active:d(S.href)}]),to:S.href,onClick:y[2]||(y[2]=F=>p.value=!1)},{default:hn(()=>[Ss(G(S.text),1)]),_:2},1032,["to","class"])]))),128))])],2),v("div",yf,[v("button",{class:"icon-btn",onClick:y[3]||(y[3]=S=>u("open-search")),onFocus:r,"aria-label":ls(t)("search")},[...y[11]||(y[11]=[v("svg",{class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[v("circle",{cx:"11",cy:"11",r:"8"}),v("line",{x1:"21",y1:"21",x2:"16.65",y2:"16.65"})],-1)])],40,vf),v("button",{class:"icon-btn",onClick:j,onFocus:r,"aria-label":ls(t)("theme")},[o.value?(E(),L("svg",xf,[...y[12]||(y[12]=[Ce('<circle cx="12" cy="12" r="5" data-v-b2e68a7c></circle><line x1="12" y1="1" x2="12" y2="3" data-v-b2e68a7c></line><line x1="12" y1="21" x2="12" y2="23" data-v-b2e68a7c></line><line x1="4.22" y1="4.22" x2="5.64" y2="5.64" data-v-b2e68a7c></line><line x1="18.36" y1="18.36" x2="19.78" y2="19.78" data-v-b2e68a7c></line><line x1="1" y1="12" x2="3" y2="12" data-v-b2e68a7c></line><line x1="21" y1="12" x2="23" y2="12" data-v-b2e68a7c></line><line x1="4.22" y1="19.78" x2="5.64" y2="18.36" data-v-b2e68a7c></line><line x1="18.36" y1="5.64" x2="19.78" y2="4.22" data-v-b2e68a7c></line>',9)])])):(E(),L("svg",Cf,[...y[13]||(y[13]=[v("path",{d:"M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"},null,-1)])]))],40,wf),v("button",{class:"icon-btn",onClick:m,onFocus:r,"aria-label":ls(t)("language")},[...y[14]||(y[14]=[Ce('<svg class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.75" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" data-v-b2e68a7c><path d="m5 8 6 6" data-v-b2e68a7c></path><path d="m4 14 6-6 2-3" data-v-b2e68a7c></path><path d="M2 5h12" data-v-b2e68a7c></path><path d="M7 2h1" data-v-b2e68a7c></path><path d="m22 22-5-10-5 10" data-v-b2e68a7c></path><path d="M14 18h6" data-v-b2e68a7c></path></svg>',1)])],40,kf)])])])])]),c.value?(E(),L("div",{key:0,class:"mobile-offcanvas d-lg-none",onClick:un(f,["self"])},[v("div",Sf,[v("div",Pf,[v("div",Af,G(ls(t)("menu")),1),v("div",Tf,[v("ul",Rf,[(E(!0),L(js,null,Ts(i.value,S=>(E(),L("li",{key:S.text,class:"my-1"},[bs(M,{to:S.href,class:Ns(["offcanvas-link d-flex align-items-center",{active:d(S.href)}]),onClick:f},{default:hn(()=>[v("span",null,G(S.text),1),y[15]||(y[15]=v("i",{class:"fas fa-chevron-right offcanvas-link__chevron"},null,-1))]),_:2},1032,["to","class"])]))),128))])])]),h.value?(E(),L("div",Ef,[v("div",Mf,G(ls(t)("tableOfContents")),1),v("div",{class:"offcanvas-tree offcanvas-card",onClick:k},[bs(qr)])])):rs("",!0)]),v("div",{class:"offcanvas-backdrop",onClick:f})])):rs("",!0)],64)}}}),Lf=mn(Df,[["__scopeId","data-v-b2e68a7c"]]),Of={class:"site-footer"},Ff={class:"container"},$f={class:"site-footer__inner"},If={class:"footer-copy"},Nf={class:"footer-designed"},Bf={class:"footer-designed__text"},zf={class:"footer-designed__name"},qf={class:"footer-designed__icons"},Hf=["href"],Vf=["href"],Gf=zs({__name:"AppFooter",setup(s){const{t:n}=Xs(),a=new Date().getFullYear(),e=en.startYear&&en.startYear<a?`${en.startYear} - ${a}`:`${a}`;return(t,l)=>(E(),L("footer",Of,[v("div",Ff,[v("div",$f,[v("p",If,"© "+G(ls(e))+" "+G(ls(en).author),1),v("div",Nf,[v("p",Bf,[Ss(G(ls(n)("designedByPrefix")),1),v("strong",zf,G(ls(en).author),1),Ss(G(ls(n)("designedBySuffix")),1)]),v("div",qf,[v("a",{href:ls(en).github,target:"_blank",rel:"noopener noreferrer",class:"footer-link","aria-label":"GitHub"},[...l[0]||(l[0]=[v("svg",{class:"footer-icon",width:"16",height:"16",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[v("path",{d:"M9 19c-5 1.5-5-2.5-7-3m14 6v-3.87a3.37 3.37 0 0 0-.94-2.61c3.14-.35 6.44-1.54 6.44-7A5.44 5.44 0 0 0 20 4.77 5.07 5.07 0 0 0 19.91 1S18.73.65 16 2.48a13.38 13.38 0 0 0-7 0C6.27.65 5.09 1 5.09 1A5.07 5.07 0 0 0 5 4.77a5.44 5.44 0 0 0-1.5 3.78c0 5.42 3.3 6.61 6.44 7A3.37 3.37 0 0 0 9 18.13V22"})],-1)])],8,Hf),v("a",{href:`mailto:${ls(en).email}`,class:"footer-link","aria-label":"Email"},[...l[1]||(l[1]=[v("svg",{class:"footer-icon",width:"16",height:"16",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[v("path",{d:"M4 4h16c1.1 0 2 .9 2 2v12c0 1.1-.9 2-2 2H4c-1.1 0-2-.9-2-2V6c0-1.1.9-2 2-2z"}),v("polyline",{points:"22,6 12,13 2,6"})],-1)])],8,Vf)])])])])]))}}),Uf=mn(Gf,[["__scopeId","data-v-f017e1f3"]]),Wf=["aria-label"],Ap="btt",Kf=zs({__name:"BackToTop",setup(s){const{t:n}=Xs(),a=ns(!1),e=ns(null),t=ns(!1),l=ns(!1),o=ns(0),p=ns(0),c=ns((typeof window<"u"?window.innerHeight:1024)-100),u=ns(!1),r=ns(null),i=ns(null);function d(){return{gap:48,minTop:68,maxTop:window.innerHeight-40-20}}function h(w){const{minTop:y,maxTop:M}=d();return Math.max(y,Math.min(M,w))}function j(w){e.value=w,!a.value&&(a.value=!0,requestAnimationFrame(()=>{window.dispatchEvent(new CustomEvent("floating-buttons-base-top",{detail:{baseTop:e.value,source:Ap}})),a.value=!1}))}function m(w){const y=w.detail,M=y?.baseTop;y?.source!==Ap&&typeof M=="number"&&(c.value=h(M))}function P(){if(r.value)return;const w=document.createElement("div");w.style.cssText="position:absolute;left:0;top:0;width:1px;height:181px;pointer-events:none;visibility:hidden",document.body.appendChild(w),i.value=w,r.value=new IntersectionObserver(([y])=>{t.value=!y.isIntersecting,t.value&&j(c.value)},{threshold:0}),r.value.observe(w)}function C(){const w=typeof window.matchMedia=="function"&&window.matchMedia("(prefers-reduced-motion: reduce)").matches;window.scrollTo({top:0,behavior:w?"auto":"smooth"})}function b(){l.value||C()}function x(w){w.preventDefault(),l.value=!0,u.value=!1,o.value=w.touches[0].clientY,p.value=c.value}function f(w){u.value=!0;const M=w.touches[0].clientY-o.value;c.value=h(p.value+M),j(c.value),w.preventDefault()}function k(w){w.preventDefault(),l.value=!1,u.value||C()}return dn(()=>{P(),window.addEventListener("floating-buttons-base-top",m),c.value=window.innerHeight-100,j(c.value)}),Tn(()=>{r.value&&(r.value.disconnect(),r.value=null),i.value&&(i.value.remove(),i.value=null),window.removeEventListener("floating-buttons-base-top",m)}),(w,y)=>Aa((E(),L("button",{class:"back-to-top d-flex align-items-center justify-content-center",onClick:b,"aria-label":ls(n)("backToTop"),onTouchstart:un(x,["prevent","stop"]),onTouchmove:un(f,["prevent","stop"]),onTouchend:un(k,["prevent","stop"]),style:oa({top:c.value+"px"})},[...y[0]||(y[0]=[v("i",{class:"fas fa-arrow-up"},null,-1)])],44,Wf)),[[sr,t.value]])}}),Yf=mn(Kf,[["__scopeId","data-v-7a6aaf37"]]),Xf={class:"page-wrap"},Jf={class:"main-content"},Qf=zs({__name:"App",setup(s){const n=ns(!1),a=iu(()=>Mp(()=>import("./SearchModal-BJTHk7sV.js"),__vite__mapDeps([0,1]))),e=Ma(),t=$l(),{locale:l}=Xs(),o=ss(()=>{const i=e.path.replace(/^\/(zh|en)(?=\/|$)/,"");return i==="/"||i===""?"":i}),p=ss(()=>l.value==="en-US"?"/en":"/zh");lt({htmlAttrs:{lang:()=>l.value==="en-US"?"en":"zh-CN"},link:ss(()=>{const i=o.value.includes("/article/"),d=i?o.value.replace(/-en$/,""):o.value,h=i&&!o.value.endsWith("-en")?`${o.value}-en`:o.value;return[{rel:"canonical",href:`${en.url}${p.value}${o.value}`},{rel:"alternate",hreflang:"zh-CN",href:`${en.url}/zh${d}`},{rel:"alternate",hreflang:"en",href:`${en.url}/en${h}`},{rel:"alternate",hreflang:"x-default",href:`${en.url}/zh${d}`}]}),meta:ss(()=>[{property:"og:url",content:`${en.url}${p.value}${o.value}`}])});const c=i=>{const d=i.match(/^\/(zh|en)(\/|$)/);d&&t.setLocale(Mr(d[1]))},u=i=>i instanceof HTMLElement?i.tagName==="INPUT"||i.tagName==="TEXTAREA"||i.isContentEditable:!1,r=i=>{i.key==="/"&&!u(i.target)&&(i.preventDefault(),n.value=!0)};return Bs(()=>e.fullPath,c),dn(()=>{c(e.fullPath),window.addEventListener("keydown",r)}),Tn(()=>{window.removeEventListener("keydown",r)}),(i,d)=>{const h=de("router-view");return E(),L("div",Xf,[bs(Lf,{onOpenSearch:d[0]||(d[0]=j=>n.value=!0)}),v("main",Jf,[bs(h,null,{default:hn(({Component:j})=>[bs(gh,{name:"view",mode:"out-in"},{default:hn(()=>[(E(),ma(yu(j)))]),_:2},1024)]),_:1})]),bs(Uf),bs(Yf),n.value?(E(),ma(ls(a),{key:0,onClose:d[1]||(d[1]=j=>n.value=!1)})):rs("",!0)])}}}),Zf=mn(Qf,[["__scopeId","data-v-73d8966b"]]);function Hr(s){return Math.max(1,Math.round(Number(s)/300))}const sb={class:"post-item__index num","aria-hidden":"true"},nb={class:"post-item__meta"},ab={class:"post-item__cat"},eb={class:"post-item__date"},tb={class:"post-item__words"},lb={class:"post-item__reading"},ob={class:"post-item__title"},pb={class:"post-item__preview"},cb={class:"post-item__tags"},rb=["onClick"],ib=["aria-label"],ub={class:"pagination__side"},hb=["aria-label"],db={class:"pagination__pages"},mb={key:1,class:"page-ellipsis"},jb=["onClick"],gb={key:2,class:"page-ellipsis"},fb={class:"pagination__side pagination__side--right"},bb=["aria-label"],_b=zs({__name:"PostList",props:{docs:{},perPage:{default:6}},setup(s){const n=s,{t:a,locale:e}=Xs(),t=Ma(),l=fe(),o=ns(1),p=ns(5),c=ns([]),u=ns([]),r=ss(()=>a("pagination")),i=ss(()=>Math.max(1,Math.ceil(n.docs.length/n.perPage))),d=ss(()=>{const O=(o.value-1)*n.perPage,V=O+n.perPage;return n.docs.slice(O,V)}),h=ss(()=>{const O=[],V=i.value,us=o.value,os=p.value;if(V<=os){for(let hs=1;hs<=V;hs++)O.push(hs);return O}O.push(us);let Y=1;for(;O.length<os&&(us-Y>=1&&!O.includes(us-Y)&&O.push(us-Y),!(O.length>=os));)us+Y<=V&&!O.includes(us+Y)&&O.push(us+Y),Y++;return O.includes(1)||(O.pop(),O.push(1)),O.length<os&&!O.includes(V)?O.push(V):O.length>=os&&!O.includes(V)&&(O.pop(),O.push(V)),O.sort((hs,Fs)=>hs-Fs)}),j=ss(()=>h.value.filter(O=>O!==1&&O!==i.value)),m=ss(()=>h.value.includes(1)),P=ss(()=>h.value.includes(i.value)),C=ss(()=>i.value>p.value&&h.value[1]>2),b=ss(()=>i.value>p.value&&h.value[h.value.length-2]<i.value-1),x=ss(()=>{const O={};try{(u.value||[]).forEach(V=>{(V.items||[]).forEach(us=>{(us.categories||[]).forEach(os=>{os&&os.key&&os.title&&(O[os.key]=os.title)})})})}catch(V){console.error(V)}return O});function f(){try{c.value=Hg()||[],u.value=ut()||[]}catch(O){console.error("Failed to load data:",O),c.value=[],u.value=[]}}function k(O){let V="";const us=Array.isArray(c.value)?c.value.find(os=>os.title===O.title):null;if(us&&us.relativePath)V=`notes/${us.relativePath}.md`;else{const os=e.value==="en-US",Y=O.category?.[1]||"notes",hs=O.title.toLowerCase().replace(/[^a-z0-9]/g,"-");V=`${Y}/${hs}.md`,os&&(V=V.replace(".md","-en.md"))}return $n(`/article/${V.replace(/\.md$/,"")}`)}function w(O){if(O>=1&&O<=i.value){o.value=O;const V={...t.query,page:String(O)};l.push({path:t.path,query:V}).catch(()=>{}),Ks(()=>window.scrollTo({top:0,behavior:"smooth"}))}}function y(){if(o.value>1){const O=o.value-1;o.value=O;const V={...t.query,page:String(O)};l.push({path:t.path,query:V}).catch(()=>{}),Ks(()=>window.scrollTo({top:0,behavior:"smooth"}))}}function M(){if(o.value<i.value){const O=o.value+1;o.value=O;const V={...t.query,page:String(O)};l.push({path:t.path,query:V}).catch(()=>{}),Ks(()=>window.scrollTo({top:0,behavior:"smooth"}))}}function S(O){if(!O)return;const V={...t.query,tag:O,page:"1"},us=e.value==="en-US"?"/en":"/zh";l.push({path:`${us}/`,query:V}).catch(()=>{}),Ks(()=>window.scrollTo({top:0,behavior:"smooth"}))}function F(){p.value=window.innerWidth<480?3:5}function Z(O){if(!Array.isArray(O)||O.length===0)return"";const V=O[0]||"",us=O[1];if(!us)return V;const os=x.value[us]||us;return`${V} / ${os}`}function B(O){return a("postReadingTime",{minutes:Hr(O)})}return f(),Bs(()=>n.docs,()=>{o.value=1}),Bs(()=>t.query.page,O=>{const V=parseInt(String(O)),us=Number.isFinite(V)&&V>=1?Math.min(V,i.value):1;us!==o.value&&(o.value=us,Ks(()=>window.scrollTo({top:0,behavior:"smooth"})))}),Bs(e,()=>f()),dn(()=>{const O=parseInt(String(t.query.page));o.value=Number.isFinite(O)&&O>=1?Math.min(O,i.value):1,F(),window.addEventListener("resize",F)}),Tn(()=>{window.removeEventListener("resize",F)}),(O,V)=>{const us=de("router-link"),os=Cl("reveal");return E(),L("div",null,[(E(!0),L(js,null,Ts(d.value,(Y,hs)=>Aa((E(),L("article",{key:Y.id,class:"post-item",style:oa({"--reveal-delay":Math.min(hs,5)*40+"ms"})},[v("span",sb,G(String((o.value-1)*n.perPage+hs+1).padStart(2,"0")),1),bs(us,{to:k(Y),class:"post-item__main"},{default:hn(()=>[v("div",nb,[v("span",ab,G(Z(Y.category)),1),V[5]||(V[5]=v("span",{class:"divider-v"},null,-1)),v("span",eb,[V[2]||(V[2]=v("i",{class:"fas fa-calendar-alt"},null,-1)),Ss(G(Y.date),1)]),V[6]||(V[6]=v("span",{class:"divider-v"},null,-1)),v("span",tb,[V[3]||(V[3]=v("i",{class:"fas fa-file-lines"},null,-1)),Ss(G(Y.wordCount)+" "+G(ls(a)("wordUnit")),1)]),V[7]||(V[7]=v("span",{class:"divider-v"},null,-1)),v("span",lb,[V[4]||(V[4]=v("i",{class:"fas fa-clock"},null,-1)),Ss(G(B(Y.wordCount)),1)])]),v("h3",ob,G(Y.title),1),v("p",pb,G(Y.preview),1)]),_:2},1032,["to"]),v("div",cb,[(E(!0),L(js,null,Ts(Y.tags,Fs=>(E(),L("span",{key:Fs,class:"post-item__tag",onClick:$s=>S(Fs)},G(Fs),9,rb))),128))]),bs(us,{to:k(Y),class:"post-item__arrow","aria-label":Y.title},{default:hn(()=>[...V[8]||(V[8]=[v("i",{class:"fas fa-arrow-right"},null,-1)])]),_:1},8,["to","aria-label"])],4)),[[os]])),128)),i.value>1?(E(),L("nav",{key:0,class:"pagination","aria-label":r.value},[v("div",ub,[o.value>1?(E(),L("button",{key:0,class:"page-btn page-btn--nav",onClick:y,"aria-label":ls(a)("prevPage")},[...V[9]||(V[9]=[v("i",{class:"fas fa-chevron-left"},null,-1)])],8,hb)):rs("",!0)]),v("div",db,[m.value?(E(),L("button",{key:0,class:Ns(["page-btn",{"page-btn--active":o.value===1}]),onClick:V[0]||(V[0]=Y=>w(1))},"1",2)):rs("",!0),C.value?(E(),L("span",mb,"...")):rs("",!0),(E(!0),L(js,null,Ts(j.value,Y=>(E(),L("button",{key:Y,class:Ns(["page-btn",{"page-btn--active":o.value===Y}]),onClick:hs=>w(Y)},G(Y),11,jb))),128)),b.value?(E(),L("span",gb,"...")):rs("",!0),P.value&&i.value>1?(E(),L("button",{key:3,class:Ns(["page-btn",{"page-btn--active":o.value===i.value}]),onClick:V[1]||(V[1]=Y=>w(i.value))},G(i.value),3)):rs("",!0)]),v("div",fb,[o.value<i.value?(E(),L("button",{key:0,class:"page-btn page-btn--nav",onClick:M,"aria-label":ls(a)("nextPage")},[...V[10]||(V[10]=[v("i",{class:"fas fa-chevron-right"},null,-1)])],8,bb)):rs("",!0)])],8,ib)):rs("",!0)])}}}),yb=mn(_b,[["__scopeId","data-v-1043b469"]]),vb={class:"page-section home-section"},wb={class:"hero"},xb={class:"hero__main"},Cb={class:"hero__greeting"},kb={class:"hero__greeting-mark"},Sb={class:"hero__name"},Pb={class:"hero__bio"},Ab={class:"hero__stats"},Tb={class:"hero__stat"},Rb={class:"hero__stat-num num"},Eb={class:"hero__stat-label"},Mb={class:"hero__stat"},Db={class:"hero__stat-num num"},Lb={class:"hero__stat-label"},Ob={class:"hero__stat"},Fb={class:"hero__stat-num num"},$b={class:"hero__stat-label"},Ib={key:0,class:"tags-row"},Nb=["onClick"],Bb={class:"posts-header"},zb={key:0,class:"posts-header__tag-title"},qb={class:"posts-header__actions"},Hb=["aria-label"],Vb=["aria-label"],Gb=zs({name:"HomeView",__name:"Home",setup(s){lt({title:"zorrooz’s blog - Home"});const{t:n,locale:a}=Xs(),e=Ma(),t=fe(),l=ns([]),o=ns([]),p=en.author,c=ss(()=>e.query.tag||""),u=ss(()=>typeof e.query.from=="string"?e.query.from:""),r=ss(()=>{const f=c.value;return f?l.value.filter(k=>Array.isArray(k.tags)&&k.tags.includes(f)):l.value}),i=ss(()=>{const f=new Map;l.value.forEach(S=>(S.tags||[]).forEach(F=>f.set(F,(f.get(F)||0)+1)));const k=Array.from(f.entries()).map(([S,F])=>({name:S,count:F})).sort((S,F)=>S.name.localeCompare(F.name)),w=k.length,y=Math.ceil(w/3),M=new Map([...k].sort((S,F)=>F.count-S.count).map((S,F)=>[S.name,F]));return k.map(S=>{const F=M.get(S.name)??0;return{...S,level:F<y?"lg":F<y*2?"md":"sm"}})}),d=ss(()=>Array.isArray(l.value)?l.value.length:0),h=ss(()=>Array.isArray(o.value)?o.value.length:0),j=ss(()=>Array.isArray(l.value)?l.value.reduce((f,k)=>{const w=typeof k?.wordCount=="number"?k.wordCount:0;return f+(Number.isFinite(w)?w:0)},0):0),m=ss(()=>{const f=j.value;return f>=1e6?(f/1e6).toFixed(f%1e6?1:0)+"M":f>=1e3?(f/1e3).toFixed(f%1e3?1:0)+"K":String(f)});function P(){try{l.value=qg()||[],o.value=Vg()||[]}catch(f){console.error("Failed to load post data:",f),l.value=[],o.value=[]}}function C(){const f={...e.query};delete f.tag,delete f.from,f.page="1",t.push({path:e.path,query:f}).catch(()=>{}),Ks(()=>{window.scrollTo({top:0,behavior:"smooth"}),P()})}function b(){u.value&&t.push(u.value).catch(()=>{})}function x(f){if(!f)return;const k={...e.query,tag:f,page:"1"},w=a.value==="en-US"?"/en":"/zh";t.push({path:`${w}/`,query:k}).catch(()=>{}),Ks(()=>window.scrollTo({top:0,behavior:"smooth"}))}return P(),Bs(a,(f,k)=>{f!==k&&(c.value?C():P())}),dn(()=>{P()}),(f,k)=>{const w=Cl("reveal");return E(),L("div",vb,[v("header",wb,[v("div",xb,[v("span",Cb,[v("span",kb,G(ls(n)("greetingPrefix")),1),Ss(G(ls(n)("greeting")),1)]),v("h1",Sb,G(ls(p)),1),v("p",Pb,G(ls(n)("developer")),1)]),v("div",Ab,[v("div",Tb,[v("span",Rb,G(d.value),1),v("span",Eb,G(ls(n)("articles")),1)]),v("div",Mb,[v("span",Db,G(h.value),1),v("span",Lb,G(ls(n)("tags")),1)]),v("div",Ob,[v("span",Fb,G(m.value),1),v("span",$b,G(ls(n)("words")),1)])])]),c.value?rs("",!0):Aa((E(),L("div",Ib,[(E(!0),L(js,null,Ts(i.value,y=>(E(),L("span",{key:y.name,class:Ns(["tag",`tag--${y.level}`]),onClick:M=>x(y.name)},G(y.name),11,Nb))),128))])),[[w]]),Aa((E(),L("div",Bb,[v("h2",{class:Ns(["posts-header__title",{"posts-header__title--tag":c.value}])},[c.value?(E(),L("span",zb,"# "+G(c.value),1)):(E(),L(js,{key:1},[Ss(G(ls(n)("recentPosts")),1)],64))],2),v("div",qb,[u.value?(E(),L("button",{key:0,class:"chip-close",onClick:b,"aria-label":ls(n)("backToArticle")},[k[0]||(k[0]=v("i",{class:"fas fa-arrow-left"},null,-1)),Ss(G(ls(n)("backToArticle")),1)],8,Hb)):rs("",!0),c.value?(E(),L("button",{key:1,class:"chip-close",onClick:C,"aria-label":ls(n)("close")},[k[1]||(k[1]=v("i",{class:"fas fa-times"},null,-1)),Ss(G(ls(n)("clearFilter")),1)],8,Vb)):rs("",!0)])])),[[w]]),bs(yb,{docs:r.value,perPage:5},null,8,["docs"])])}}}),Ub=mn(Gb,[["__scopeId","data-v-f7cc79bf"]]),Wb={class:"page-section category-view"},Kb={class:"category-head"},Yb={class:"article-title"},Xb={class:"cat-section__header"},Jb={class:"cat-section__title"},Qb={class:"cat-section__count"},Zb={class:"cat-grid"},s_={class:"cat-card__head"},n_={class:"cat-card__name"},a_={class:"cat-card__ext-links"},e_=["href"],t_=["href"],l_={class:"cat-card__desc"},o_={key:0,class:"cat-card__stats"},p_={key:0,class:"cat-stat"},c_={key:1,class:"cat-stat"},r_={key:2,class:"cat-stat"},i_={key:1,class:"cat-card__stats"},u_={key:0,class:"cat-stat"},h_={key:1,class:"cat-stat"},d_={key:2,class:"cat-stat"},m_={key:2,class:"cat-card__tags"},j_={class:"cat-card__links"},g_=["onClick"],f_=zs({name:"CategoryView",__name:"Category",setup(s){lt({title:"zorrooz’s blog - Categories"});const{t:n,locale:a}=Xs(),e=fe(),t=ns([]),l=ss(()=>n("categories")),o=ss(()=>n("notes")),p=ss(()=>n("projects")),c=ss(()=>n("topics")),u=ss(()=>n("seeMore"));function r(){try{const b=ut();t.value=Array.isArray(b)?b:b?.categoryList||[]}catch(b){console.error("Failed to load category data:",b),t.value=[]}}function i(b){return b===o.value?"fa-book-open":b===p.value?"fa-folder-open":b===c.value?"fa-flask":"fa-folder"}function d(b){const x=b?.items||[];if(b.title===o.value){const f=x.reduce((k,w)=>k+(w.stats?.postsCount||0),0);return n("countPosts",{count:f})}return b.title===p.value?n("countProjects",{count:x.length}):n("countTopics",{count:x.length})}function h(b){return b&&b.stats&&b.stats.latestDate||""}function j(b){return typeof b!="string"||!b.trim()?"":/^https?:\/\//i.test(b)?b:"https://"+b.replace(/^\/+/,"")}function m(b){return typeof b!="string"||!b.trim()?"":/^https?:\/\//i.test(b)?b:/^10\.\d{4,9}\//.test(b)?"https://doi.org/"+b:"https://"+b.replace(/^\/+/,"")}function P(b){return!!(b?.url||b?.github||b?.doi)}function C(b){const x=b?.root;if(x){e.push($n(x)).catch(y=>{y.name!=="NavigationDuplicated"&&!y.toString().includes("Navigation cancelled")&&console.error("Navigation error:",y)});return}const f=j(b?.url);if(f){window.open(f,"_blank","noopener,noreferrer");return}const k=j(b?.github);if(k){window.open(k,"_blank","noopener,noreferrer");return}const w=m(b?.doi);w&&window.open(w,"_blank","noopener,noreferrer")}return r(),Bs(a,()=>{r()}),(b,x)=>{const f=Cl("reveal");return E(),L("div",Wb,[v("div",Kb,[v("h1",Yb,G(l.value),1)]),(E(!0),L(js,null,Ts(t.value,(k,w)=>(E(),L("div",{key:w,class:"cat-section"},[v("div",Xb,[v("h2",Jb,[v("i",{class:Ns([["fas",i(k.title)],"cat-section__icon"]),"aria-hidden":"true"},null,2),Ss(" "+G(k.title),1)]),v("span",Qb,G(d(k)),1)]),v("div",Zb,[(E(!0),L(js,null,Ts(k.items,(y,M)=>Aa((E(),L("div",{key:M,class:"cat-card",style:oa({"--reveal-delay":Math.min(Number(M),5)*40+"ms"})},[v("div",s_,[v("div",n_,G(y.title),1),v("div",a_,[y.github?(E(),L("a",{key:0,href:j(y.github),target:"_blank",rel:"noopener noreferrer",class:"cat-card__ext-link","aria-label":"GitHub"},[...x[0]||(x[0]=[v("i",{class:"fab fa-github"},null,-1)])],8,e_)):rs("",!0),y.doi?(E(),L("a",{key:1,href:m(y.doi),target:"_blank",rel:"noopener noreferrer",class:"cat-card__ext-link","aria-label":"DOI"},[...x[1]||(x[1]=[v("i",{class:"fas fa-link"},null,-1)])],8,t_)):rs("",!0)])]),v("p",l_,G(y.desc),1),k.title===o.value?(E(),L("div",o_,[y.stats?.postsCount?(E(),L("span",p_,[x[2]||(x[2]=v("i",{class:"fas fa-file-lines"},null,-1)),Ss(G(ls(n)("countPosts",{count:y.stats.postsCount})),1)])):rs("",!0),y.stats?.totalWords?(E(),L("span",c_,[x[3]||(x[3]=v("i",{class:"fas fa-font"},null,-1)),Ss(G(ls(n)("countWords",{count:y.stats.totalWords})),1)])):rs("",!0),h(y)?(E(),L("span",r_,[x[4]||(x[4]=v("i",{class:"fas fa-clock"},null,-1)),Ss(G(h(y)),1)])):rs("",!0)])):(E(),L("div",i_,[y.language?(E(),L("span",u_,[x[5]||(x[5]=v("i",{class:"fas fa-code"},null,-1)),Ss(G(y.language),1)])):rs("",!0),y.year?(E(),L("span",h_,[x[6]||(x[6]=v("i",{class:"fas fa-calendar"},null,-1)),Ss(G(y.year),1)])):rs("",!0),y.license?(E(),L("span",d_,[x[7]||(x[7]=v("i",{class:"fas fa-scale-balanced"},null,-1)),Ss(G(y.license),1)])):rs("",!0)])),Array.isArray(y.tags)&&y.tags.length?(E(),L("div",m_,[(E(!0),L(js,null,Ts(y.tags,(S,F)=>(E(),L("span",{key:F,class:"cat-card__tag"},G(S),1))),128))])):rs("",!0),v("div",j_,[P(y)||y.root?(E(),L("a",{key:0,class:"cat-card__link",onClick:un(S=>C(y),["prevent"])},[Ss(G(u.value),1),x[8]||(x[8]=v("i",{class:"fas fa-arrow-right"},null,-1))],8,g_)):rs("",!0)])],4)),[[f]])),128))])]))),128))])}}}),b_=mn(f_,[["__scopeId","data-v-fd52f9cf"]]),__={class:"page-section resource-view"},y_={class:"resource-head"},v_={class:"article-title"},w_={class:"resource-subtitle"},x_={class:"res-layout"},C_={class:"res-sidebar"},k_={class:"res-group__label"},S_={class:"res-group__count num"},P_={class:"res-group__items"},A_=["onClick"],T_={class:"res-main"},R_={class:"res-groups"},E_={class:"res-group-block__title"},M_={class:"res-grid"},D_=["href"],L_={class:"res-card__head"},O_={class:"res-card__name"},F_={class:"res-card__ext-links"},$_={key:0,class:"res-card__ext-link","aria-label":"DOI"},I_={class:"res-card__desc"},N_={class:"res-card__footer"},B_={class:"res-card__url"},z_=zs({name:"ResourceView",__name:"Resource",setup(s){const{t:n,locale:a}=Xs(),e=ns([]),t=ns(null),l=ss(()=>n("resources")),o=ss(()=>n("resourceSubtitle")),p=ss(()=>t.value?Array.isArray(t.value.children)&&t.value.children.length?t.value.children.filter(h=>Array.isArray(h.items)&&h.items.length):Array.isArray(t.value.items)?[{title:t.value.title,items:t.value.items}]:[]:[]);function c(h){t.value=h}function u(h){return t.value===h}function r(h){return h?h.replace(/^https?:\/\//,"").replace(/\/$/,""):""}function i(h){return!!h&&(h.includes("doi.org")||/^10\.\d{4,9}\//.test(h))}function d(){try{e.value=Ug()||[],t.value=e.value[0]?.children?.[0]||null}catch(h){console.error("Failed to load resources data:",h),e.value=[],t.value=null}}return d(),Bs(a,()=>{d()}),(h,j)=>(E(),L("div",__,[v("header",y_,[v("h1",v_,G(l.value),1),v("p",w_,[j[0]||(j[0]=v("i",{class:"fas fa-circle-info resource-head__icon"},null,-1)),Ss(G(o.value),1)])]),v("div",x_,[v("aside",C_,[(E(!0),L(js,null,Ts(e.value,m=>(E(),L("div",{key:m.title,class:"res-group"},[v("div",k_,[v("span",null,G(m.title),1),v("span",S_,G(m.children?.length||0),1)]),v("div",P_,[(E(!0),L(js,null,Ts(m.children,P=>(E(),L("button",{key:P.title,class:Ns(["res-item",{"res-item--active":u(P)}]),onClick:C=>c(P)},G(P.title),11,A_))),128))])]))),128))]),j[3]||(j[3]=v("div",{class:"res-divider","aria-hidden":"true"},null,-1)),v("main",T_,[v("div",R_,[(E(!0),L(js,null,Ts(p.value,m=>(E(),L("section",{key:m.title,class:"res-group-block"},[v("h3",E_,G(m.title),1),v("div",M_,[(E(!0),L(js,null,Ts(m.items,P=>(E(),L("a",{key:P.name,href:P.url,target:"_blank",rel:"noopener noreferrer",class:"res-card"},[v("div",L_,[v("span",O_,G(P.name),1),v("span",F_,[i(P.url)?(E(),L("span",$_,[...j[1]||(j[1]=[v("i",{class:"fas fa-link"},null,-1)])])):rs("",!0)])]),v("p",I_,G(P.desc),1),v("div",N_,[v("span",B_,G(r(P.url)),1),j[2]||(j[2]=v("i",{class:"fas fa-arrow-up-right-from-square res-card__arrow"},null,-1))])],8,D_))),128))])]))),128))])])])]))}}),q_=mn(z_,[["__scopeId","data-v-fdc07a4f"]]),H_="/assets/avatar-DQvqWlfS.png",V_={class:"page-section about-view"},G_={class:"about-head"},U_={class:"about-head__identity"},W_={class:"about-head__avatar"},K_=["src"],Y_={key:1,class:"about-head__initial"},X_={class:"about-head__names"},J_={class:"about-head__name"},Q_={class:"about-head__role"},Z_={key:0,class:"about-intro"},sy={class:"about-body"},ny={key:0,class:"about-section"},ay={class:"about-section__title"},ey={class:"timeline"},ty={class:"tl-year num"},ly={class:"tl-content"},oy={class:"tl-title"},py={key:0,class:"tl-desc"},cy={key:1,class:"about-grid"},ry={class:"about-cell__title"},iy={key:0,class:"about-cell__item-name"},uy={key:1,class:"about-cell__item-desc"},hy={class:"about-foot"},dy={class:"about-foot__contacts"},my=["href"],jy={class:"about-foot__thanks"},gy=zs({name:"AboutView",__name:"About",setup(s){const{t:n,locale:a}=Xs(),e=ns({}),t=Object.assign({"/src/assets/avatar.png":H_}),l=ss(()=>Object.values(t)[0]||""),o=ss(()=>n("thanks")),p=en.author,c=p.trim().charAt(0).toUpperCase(),u=ss(()=>e.value?.introduction||""),r=ss(()=>e.value?.experience||[]),i=ss(()=>e.value?.section||[]);function d(){try{e.value=Gg()||{}}catch(h){console.error("Failed to load about data:",h),e.value={introduction:"",experience:[],section:[],contacts:[]}}}return d(),Bs(a,()=>{d()}),(h,j)=>(E(),L("div",V_,[v("header",G_,[v("div",U_,[v("div",W_,[l.value?(E(),L("img",{key:0,src:l.value,alt:"avatar",class:"about-head__avatar-img",draggable:"false"},null,8,K_)):(E(),L("span",Y_,G(ls(c)),1))]),v("div",X_,[v("h1",J_,G(ls(p)),1),v("p",Q_,G(ls(n)("developer")),1)])]),u.value?(E(),L("p",Z_,G(u.value),1)):rs("",!0)]),v("main",sy,[r.value.length?(E(),L("section",ny,[v("div",ay,G(ls(n)("experience")),1),v("div",ey,[(E(!0),L(js,null,Ts(r.value,(m,P)=>(E(),L("div",{key:P,class:"tl-item"},[v("div",ty,G(m.year),1),v("div",ly,[v("div",oy,G(m.title),1),m.desc?(E(),L("div",py,G(m.desc),1)):rs("",!0)])]))),128))])])):rs("",!0),i.value.length?(E(),L("div",cy,[(E(!0),L(js,null,Ts(i.value,(m,P)=>(E(),L("div",{key:P,class:"about-cell"},[v("div",ry,G(m.title),1),(E(!0),L(js,null,Ts(m.items,(C,b)=>(E(),L("div",{key:b,class:"about-cell__item"},[C.name?(E(),L("span",iy,G(C.name),1)):rs("",!0),C.desc?(E(),L("span",uy,G(C.desc),1)):rs("",!0)]))),128))]))),128))])):rs("",!0)]),v("footer",hy,[v("div",dy,[(E(!0),L(js,null,Ts(e.value.contacts,m=>(E(),L("a",{key:m.label,href:m.link,target:"_blank",rel:"noopener noreferrer",class:"about-foot__contact"},[v("i",{class:Ns(m.icon)},null,2),v("span",null,G(m.value),1)],8,my))),128))]),v("p",jy,G(o.value),1)])]))}}),fy=mn(gy,[["__scopeId","data-v-4ad5ac55"]]),by=["innerHTML"],Tp=`
  <svg width="16" height="16" viewBox="0 0 14 14" fill="currentColor">
    <path d="M3 2C2.44772 2 2 2.44772 2 3V9C2 9.55228 2.44772 10 3 10H9C9.55228 10 10 9.55228 10 9V3C10 2.44772 9.55228 2 9 2H3ZM1 3C1 1.89543 1.89543 1 3 1H9C10.1046 1 11 1.89543 11 3V9C11 10.1046 10.1046 11 9 11H3C1.89543 11 1 10.1046 1 9V3Z"/>
    <path d="M5 4C4.44772 4 4 4.44772 4 5V11C4 11.5523 4.44772 12 5 12H11C11.5523 12 12 11.5523 12 11V5C12 4.44772 11.5523 4 11 4H5Z"/>
  </svg>
`,_y=`
  <svg width="16" height="16" viewBox="0 0 14 14" fill="currentColor">
    <path d="M11.3536 3.64645C11.5488 3.84171 11.5488 4.15829 11.3536 4.35355L5.35355 10.3536C5.15829 10.5488 4.84171 10.5488 4.64645 10.3536L2.64645 8.35355C2.45118 8.15829 2.45118 7.84171 2.64645 7.64645C2.84171 7.45118 3.15829 7.45118 3.35355 7.64645L5 9.29289L10.6464 3.64645C10.8417 3.45118 11.1583 3.45118 11.3536 3.64645Z"/>
  </svg>
`,yy=zs({__name:"RenderMarkdown",props:{rawMarkdown:{default:""},articlePath:{default:""},articleTitle:{default:""}},emits:["markdown-rendered"],setup(s,{emit:n}){const{t:a}=Xs(),e=Object.assign({}),t=s,l=n,o=ns(""),p=Te("markdownContainer");async function c(b){const x=u(b,t.articlePath);o.value=x,await Ks(),l("markdown-rendered"),h(),j(),r()}function u(b,x){try{const f=x.replace(/^[./]*/,"").replace(/\.md$/,"").split("/").slice(0,-1).join("/"),k=w=>{if(/^(https?:)?\/\//i.test(w)||w.startsWith("/"))return w;const y=(f+"/"+w).split("/").filter(Z=>Z&&Z!=="."),M=[];y.forEach(Z=>Z===".."?M.pop():M.push(Z));const S=M.join("/"),F=[`../../content-src/${S}`,`../content-src/${S}`,`content-src/${S}`];for(const Z of F)if(e[Z])return e[Z];return w};return b.replace(/<img\s+([^>]*?)src=["']([^"']+)["'](.*?)>/gi,(w,y,M,S)=>`<img ${y}src="${k(M.trim())}"${S}>`)}catch(f){return console.warn("rewriteImageLinks failed",f),b}}function r(){const b=p.value;b&&(i(b),d(b))}function i(b){if(!t.articleTitle)return;const x=t.articleTitle.trim().toLowerCase();b.querySelectorAll("h1").forEach(f=>{f.textContent.trim().toLowerCase()===x?f.remove():f.replaceWith(Object.assign(document.createElement("h2"),{...Object.fromEntries(Array.from(f.attributes).map(w=>[w.name,w.value])),innerHTML:f.innerHTML}))})}function d(b){const x=(f,k)=>{const w=window.scrollY+f.getBoundingClientRect().top-88;window.scrollTo({top:Math.max(0,w),behavior:"smooth"}),setTimeout(()=>k.blur(),300)};b.querySelectorAll("h2, h3, h4, h5, h6").forEach(f=>{f.querySelector(".heading-anchor")?.remove();const k=Object.assign(document.createElement("button"),{type:"button",className:"heading-anchor",textContent:"#",ariaLabel:a("anchorHeading"),tabIndex:0,ariaHidden:"false"});k.addEventListener("click",w=>{w.stopPropagation(),x(f,k)}),k.addEventListener("keydown",w=>{(w.key==="Enter"||w.key===" ")&&(w.preventDefault(),x(f,k))}),f.appendChild(k)})}function h(){const b=p.value;b&&b.querySelectorAll("pre").forEach(x=>{if(x.querySelector(".code-block-header"))return;const f=x.querySelector("code");if(!f)return;const k=(f.className.match(/language-(\w+)/)||["","text"])[1],w=document.createElement("div");w.className="code-block-header d-flex align-items-center justify-content-between";const y=document.createElement("span");y.className="code-language",y.textContent=k;const M=document.createElement("button");M.type="button",M.className="copy-button btn-icon d-flex align-items-center justify-content-center",M.setAttribute("aria-label",a("copyCode")),M.innerHTML=Tp,M.addEventListener("click",()=>P(f.textContent??"",M)),w.append(y,M);const S=document.createElement("div");S.className="code-block-wrapper",x.parentNode?.insertBefore(S,x),S.append(w,x)})}function j(){const b=p.value;b&&b.querySelectorAll("table").forEach(x=>{if(x.closest(".table-copyable"))return;const f=document.createElement("div");f.className="table-copyable";const k=document.createElement("button");k.type="button",k.className="table-copy-btn",k.setAttribute("aria-label",a("copyTable")),k.innerHTML=Tp,k.addEventListener("click",()=>m(x,k)),f.append(k),x.parentNode?.insertBefore(f,x),f.append(x)})}function m(b,x){const k=Array.from(b.querySelectorAll("tr")).map(w=>Array.from(w.querySelectorAll("th, td")).map(y=>(y.textContent||"").trim().replace(/\s+/g," ")).join("	"));P(k.join(`
`),x)}async function P(b,x){try{await navigator.clipboard.writeText(b)}catch(f){console.error(a("copyFailed"),f);const k=document.createElement("textarea");k.value=b,document.body.appendChild(k),k.select(),document.execCommand("copy"),document.body.removeChild(k)}finally{C(x)}}function C(b){const x=b.innerHTML;b.style.color="var(--primary)",b.innerHTML=_y,setTimeout(()=>{b.innerHTML=x,b.style.color=""},1200)}return Bs(()=>t.rawMarkdown,b=>c(b),{immediate:!0}),(b,x)=>(E(),L("div",{class:"markdown-body",innerHTML:o.value,ref_key:"markdownContainer",ref:p},null,8,by))}}),vy={class:"on-this-page"},wy={class:"otp-header"},xy={class:"otp-title"},Cy={class:"otp-list"},ky=["onClick","onKeydown"],Sy={class:"otp-text"},Py={key:0,class:"otp-sublist"},Ay=["onClick","onKeydown"],Ty={class:"otp-subtext"},Ry=zs({__name:"OnThisPage",props:{containerSelector:{default:".markdown-body"},levels:{default:()=>[2,3,4,5,6]},offset:{default:8}},emits:["navigate"],setup(s,{expose:n,emit:a}){const e=s,t=a,{t:l}=Xs(),o=ns([]),p=ns(""),c=ns(null),u=ns(null),r=ns(null),i=ns(null);function d(){c.value&&(c.value.disconnect(),c.value=null),u.value&&(clearTimeout(u.value),u.value=null),r.value&&(clearInterval(r.value),r.value=null),i.value&&(i.value.disconnect(),i.value=null)}function h(){o.value=[],p.value="",d(),Ks(()=>m())}function j(){C(),b()}function m(){setTimeout(()=>j(),100);const f=()=>{const k=document.querySelector(e.containerSelector);k&&(r.value&&(clearInterval(r.value),r.value=null),c.value||(c.value=new MutationObserver(()=>{clearTimeout(u.value??void 0),u.value=window.setTimeout(()=>j(),100)}),c.value.observe(k,{childList:!0,subtree:!0,attributes:!0,attributeFilter:["id"]})),j())};f(),!r.value&&!document.querySelector(e.containerSelector)&&(r.value=window.setInterval(f,200))}function P(f){try{const k=f.cloneNode(!0);return k.querySelectorAll(".heading-anchor")?.forEach(w=>w.remove()),(k.textContent||"").replace(/\s*#\s*$/,"").trim()}catch{return(f.textContent||"").replace(/\s*#\s*$/,"").trim()}}function C(){const f=document.querySelector(e.containerSelector);if(!f){o.value=[];return}const k=e.levels.map(B=>`h${B}`).join(","),w=Array.from(f.querySelectorAll(k));if(w.length===0){o.value=[];return}w.forEach(B=>{if(!B.id){let O=B.textContent.trim().toLowerCase().replace(/[^\u4e00-\u9fa5a-zA-Z0-9\s-]/g,"").replace(/\s+/g,"-");O||(O=`section-${Math.random().toString(36).substring(2,9)}`);let V=O,us=1;for(;document.getElementById(V);)V=`${O}-${us++}`;B.id=V}});const y=new Set(e.levels),M=Math.min(...e.levels),S=M+1,F=[];let Z=null;for(const B of w){const O=parseInt(B.tagName.substring(1),10);if(!y.has(O))continue;const V={id:B.id,text:P(B),level:O,children:[]};O===M?(F.push(V),Z=V):Z&&O>=S?Z.children.push(V):F.push(V)}o.value=F}function b(){i.value&&(i.value.disconnect(),i.value=null);const f=document.querySelector(e.containerSelector);if(!f)return;const k=e.levels.map(y=>`h${y}`).join(","),w=Array.from(f.querySelectorAll(k));w.length!==0&&(i.value=new IntersectionObserver(y=>{for(const M of y)M.isIntersecting&&(p.value=M.target.id)},{rootMargin:`-${e.offset}px 0px -60% 0px`,threshold:0}),w.forEach(y=>i.value?.observe(y)))}function x(f){t("navigate",f);const k=document.getElementById(f);if(!k)return;const w=k.getBoundingClientRect().top+window.scrollY-e.offset,y=typeof window.matchMedia=="function"&&window.matchMedia("(prefers-reduced-motion: reduce)").matches,M=()=>{window.scrollTo({top:w,behavior:y?"auto":"smooth"})};try{(document.body&&document.body.style&&document.body.style.overflow)==="hidden"?setTimeout(M,80):M()}catch{M()}}return dn(()=>{C(),b(),Ks(()=>{m()})}),Tn(()=>{d()}),n({refreshToc:j,resetToc:h}),(f,k)=>(E(),L("nav",vy,[v("div",wy,[v("span",xy,G(ls(l)("tableOfContents")),1)]),v("ul",Cy,[(E(!0),L(js,null,Ts(o.value,(w,y)=>(E(),L("li",{key:y,class:Ns(["otp-item",{active:p.value===w.id}])},[v("a",{class:"otp-link",role:"button",tabindex:"0",onClick:un(M=>x(w.id),["prevent"]),onKeydown:Lo(un(M=>x(w.id),["prevent"]),["enter"])},[v("span",Sy,G(w.text),1)],40,ky),w.children&&w.children.length?(E(),L("ul",Py,[(E(!0),L(js,null,Ts(w.children,(M,S)=>(E(),L("li",{key:S,class:Ns(["otp-subitem",{active:p.value===M.id}])},[v("a",{class:"otp-sublink",role:"button",tabindex:"0",onClick:un(F=>x(M.id),["prevent"]),onKeydown:Lo(un(F=>x(M.id),["prevent"]),["enter"])},[v("span",Ty,G(M.text),1)],40,Ay)],2))),128))])):rs("",!0)],2))),128))])]))}}),Vr=mn(Ry,[["__scopeId","data-v-309aab29"]]),Ey=["aria-label"],My={class:"offcanvas-panel offcanvas-right border-start rounded-0 shadow-sm"},Dy={class:"offcanvas-section"},Ly={class:"offcanvas-card"},Rp="toc",Oy=zs({__name:"TocDrawer",setup(s){const{t:n}=Xs(),a=ns(!1),e=ns(null),t=ns(!1),l=ns(!1),o=ns(0),p=ns(0),c=ns((typeof window<"u"?window.innerHeight:1024)-160),u=ns(!1),r=ns(!1),i=ns(null);function d(){return{gap:48,minTop:20,maxTop:Math.max(0,window.innerHeight-40-20-48)}}function h(S){const{minTop:F,maxTop:Z}=d();return Math.max(F,Math.min(Z,S))}function j(S){e.value=S,!a.value&&(a.value=!0,requestAnimationFrame(()=>{window.dispatchEvent(new CustomEvent("floating-buttons-base-top",{detail:{baseTop:e.value,source:Rp}})),a.value=!1}))}function m(S){const F=S.detail,Z=F?.baseTop;if(F?.source===Rp)return;const{gap:O}=d();if(typeof Z=="number"){const V=Z-O;c.value=h(V)}}function P(){const S=window.innerWidth<992;t.value=S,c.value=h(c.value)}function C(){l.value||(y(),r.value=!0)}function b(){r.value=!1,M()}function x(S){S.preventDefault(),l.value=!0,u.value=!1,o.value=S.touches[0].clientY,p.value=c.value}function f(S){u.value=!0;const Z=S.touches[0].clientY-o.value,B=h(p.value+Z);c.value=B;const{gap:O}=d();j(B+O),S.preventDefault()}function k(S){S.preventDefault(),l.value=!1,u.value||C()}function w(){b()}function y(){try{const S=window.scrollY||window.pageYOffset||document.documentElement.scrollTop||0;i.value=S;const F=document.body;F&&(F.style.position="fixed",F.style.top=`-${S}px`,F.style.left="0",F.style.right="0",F.style.overflow="hidden"),document.documentElement&&(document.documentElement.style.overscrollBehavior="contain")}catch{const S=document.documentElement,F=document.body;S&&(S.style.overflow="hidden",S.style.overscrollBehavior="contain"),F&&(F.style.overflow="hidden",F.style.overscrollBehavior="contain")}}function M(){try{const S=document.body;S&&(S.style.position="",S.style.top="",S.style.left="",S.style.right="",S.style.overflow=""),document.documentElement&&(document.documentElement.style.overscrollBehavior=""),typeof i.value=="number"&&(window.scrollTo(0,i.value),i.value=null)}catch{const S=document.documentElement,F=document.body;S&&(S.style.overflow="",S.style.overscrollBehavior=""),F&&(F.style.overflow="",F.style.overscrollBehavior="")}}return dn(()=>{window.addEventListener("resize",P,{passive:!0}),window.addEventListener("floating-buttons-base-top",m),P(),j(c.value+48)}),Tn(()=>{window.removeEventListener("resize",P),window.removeEventListener("floating-buttons-base-top",m),M()}),(S,F)=>(E(),L(js,null,[Aa(v("button",{class:"toc-drawer-btn d-lg-none d-flex align-items-center justify-content-center",onClick:C,"aria-label":ls(n)("openToc"),onTouchstart:un(x,["prevent","stop"]),onTouchmove:un(f,["prevent","stop"]),onTouchend:un(k,["prevent","stop"]),style:oa({top:c.value+"px"})},[...F[0]||(F[0]=[v("i",{class:"fas fa-bookmark"},null,-1)])],44,Ey),[[sr,t.value]]),r.value?(E(),L("div",{key:0,class:"mobile-offcanvas d-lg-none",onClick:un(b,["self"])},[v("div",My,[v("div",Dy,[v("div",Ly,[bs(Vr,{containerSelector:".markdown-body",levels:[2,3],offset:88,onNavigate:w})])])]),v("div",{class:"offcanvas-backdrop",onClick:b})])):rs("",!0)],64))}}),Fy=mn(Oy,[["__scopeId","data-v-3f1455a9"]]),$y={class:"container view-container article-view"},Iy={class:"row py-4 px-0"},Ny={class:"col-12 col-lg-2 order-2 order-lg-1 docs-sidebar-col"},By={key:0,class:"navigation-container mb-0"},zy={class:"col-12 col-lg-8 order-1 order-lg-2 docs-main-col",ref:"mainContent"},qy={class:"article-panel"},Hy={class:"article-panel__body"},Vy={class:"article-content"},Gy={key:0,class:"article-meta"},Uy={class:"article-title"},Wy={class:"article-meta__row"},Ky={key:0,class:"meta-line"},Yy={key:1,class:"meta-line"},Xy=["aria-label"],Jy={key:0,class:"article-meta__tags"},Qy=["onClick"],Zy={key:2,class:"article-navigation"},sv={key:0,class:"article-nav-spacer","aria-hidden":"true"},nv={class:"nav-details"},av={class:"nav-label"},ev={class:"nav-title"},tv={class:"nav-details"},lv={class:"nav-label"},ov={class:"nav-title"},pv={class:"col-12 col-lg-2 order-3 docs-toc-col"},cv={key:0,class:"toc-container mt-0"},rv=zs({name:"ArticleView",head(){return{title:this.currentPost?`${this.currentPost.title} - zorrooz’s blog`:"zorrooz’s blog",meta:this.currentPost?.description?[{name:"description",content:this.currentPost.description}]:[]}},__name:"Article",props:{path:{type:[String,Array],default:""}},setup(s){const{t:n,locale:a}=Xs(),e=Ma(),t=fe(),l=ns(""),o=ns(""),p=ns([]),c=ns([]),u=ns(typeof window<"u"?window.innerWidth:1024),r=ns(0),i=ns(!1);let d=null,h=!1;const j=Te("onThisPageRef"),m=Te("leftSidebarContent"),P=Te("rightSidebarContent"),C=ss(()=>u.value>=992),b=ss(()=>!!(w.value?.path&&w.value.path.startsWith("notes/"))),x=ss(()=>n("updatedAt")),f=ss(()=>n("prevPage")),k=ss(()=>n("nextPage")),w=ss(()=>o.value?p.value.find(es=>{const ms=es.path.replace(/\.md$/,""),_s=o.value.replace(/\.md$/,"");return!!(ms===_s||_s.endsWith("-en")&&ms===_s.replace(/-en$/,"")||!_s.endsWith("-en")&&ms===_s+"-en")}):null),y=ss(()=>{if(!w.value)return[];const[es,ms]=w.value.path.replace(/\.md$/,"").split("/"),_s=[],fs=(R,H)=>{if(!H)return;const q=String(H).replace(/^\/+/,"").split("/"),X=q[0]==="article"?1:0,gs=q[X],g=q[X+1];if(gs!==es||g!==ms)return;const _=q.slice(X+2),A=`${gs}/${g}/${_.join("/")}`;_s.push({title:R,path:`${A}.md`})};if(Array.isArray(c.value)){for(const R of c.value)if(Array.isArray(R.items))for(const H of R.items)H?.name===ms&&(Array.isArray(H.articles)&&H.articles.forEach(q=>fs(q.title,q.articleUrl)),Array.isArray(H.categories)&&H.categories.forEach(q=>{Array.isArray(q.articles)&&q.articles.forEach(X=>fs(X.title,X.articleUrl))}))}return _s}),M=ss(()=>w.value?y.value.findIndex(es=>es.path.replace(/\.md$/,"")===w.value.path.replace(/\.md$/,"")):-1),S=ss(()=>{const es=M.value;return es>0?y.value[es-1]:null}),F=ss(()=>{const es=M.value,ms=y.value.length-1;return es>=0&&es<ms?y.value[es+1]:null}),Z=ss(()=>{const es=w.value?.wordCount;return typeof es=="number"?Hr(es):0});function B(es){return n("articleReadingTime",{minutes:es})}function O(es){return $n(`/article/${es.replace(/\.md$/,"")}`)}function V(es){if(!es)return;const ms=a.value==="en-US"?"/en":"/zh";t.push({path:`${ms}/`,query:{tag:es,from:e.fullPath}})}async function us(){const es=document.querySelector(".markdown-body");if(!es||!w.value)return;const ms=es.cloneNode(!0);ms.querySelectorAll(".heading-anchor, .code-block-header, .copy-button, .table-copy-btn").forEach(R=>R.remove());const _s=(ms.innerText||"").trim(),fs=`${w.value.title}

${_s}`;try{await navigator.clipboard.writeText(fs)}catch(R){console.error(n("copyFailed"),R);const H=document.createElement("textarea");H.value=fs,document.body.appendChild(H),H.select(),document.execCommand("copy"),document.body.removeChild(H)}finally{i.value=!0,d&&clearTimeout(d),d=window.setTimeout(()=>{i.value=!1},1200)}}function os(){const es=[],ms=(_s,fs,R=[],H="",q=0)=>{if(typeof fs!="string"||!fs.trim())return;const X=fs.replace(/^\/+/,"").split("/"),gs=X[0]==="article"?1:0,g=X[gs],_=X[gs+1],A=X.slice(gs+2),$=A[0]||"__root__",N=A.join("/"),I=`${g}/${_}/${N}`,K={title:_s,path:`${I}.md`,date:H||"",tags:Array.isArray(R)?R:[],preview:"",category:`${g}/${_}/${$}`,wordCount:typeof q=="number"?q:0};es.push(K)};try{const _s=ut();Array.isArray(_s)&&_s.forEach(fs=>{Array.isArray(fs.items)&&fs.items.forEach(R=>{const H=R?.stats?.latestDate||"";Array.isArray(R.articles)&&R.articles.forEach(q=>ms(q.title,q.articleUrl,q?.tags||[],H,q?.wordCount)),Array.isArray(R.categories)&&R.categories.forEach(q=>{const X=q?.stats?.latestDate||H||"";Array.isArray(q.articles)&&q.articles.forEach(gs=>ms(gs.title,gs.articleUrl,gs?.tags||[],X,gs?.wordCount))})})}),p.value=es,c.value=_s}catch{p.value=[],c.value=[]}}function Y(){u.value=window.innerWidth,Hs()}function hs(){h||(h=!0,requestAnimationFrame(()=>{h=!1;const ms=document.documentElement.scrollHeight-window.innerHeight;r.value=ms>0?Math.min(100,window.scrollY/ms*100):0}))}function Fs(es){return Array.isArray(es)?es.join("/"):typeof es=="string"&&es.length?es:""}function $s(){l.value="";try{const es=Fs(e.params.path),ms=a.value==="en-US";let _s=p.value.find(fs=>{const R=fs.path.replace(/\.md$/,"");return ms?R===es||R===es.replace(/-en$/,"")+"-en":R===es||R===es.replace(/-en$/,"")});if(_s||(_s=p.value.find(fs=>{const R=fs.path.replace(/\.md$/,""),H=es.replace(/-en$/,"");return R===H||R===H+"-en"})),!_s)throw new Error(`Article not found: ${es}`);o.value=_s.path,l.value=zg(o.value),Ks(()=>{typeof window>"u"||(window.scrollTo({top:0,behavior:"smooth"}),Hs(),j.value?.refreshToc())})}catch{l.value=`# Article Not Found

The requested article could not be loaded. Please check the URL.`,Ks(()=>{typeof window>"u"||(window.scrollTo({top:0,behavior:"smooth"}),j.value?.refreshToc())})}}function Hs(){if(typeof window>"u")return;const ms=document.querySelector("header")?.offsetHeight||60,_s=window.innerHeight,fs=Math.max(200,_s-ms-24-24);[m.value,P.value].forEach(H=>{H&&(H.style.maxHeight=`${fs}px`,H.style.overflowY="auto")})}function Vs(){Hs(),Ks(()=>j.value?.refreshToc())}return os(),$s(),Bs(a,(es,ms)=>{es!==ms&&(os(),$s())}),Bs(()=>e.params.path,(es,ms)=>{const _s=Fs(ms),fs=Fs(es);_s!==fs&&(j.value?.resetToc(),$s())}),Bs(l,()=>{Ks(()=>Hs())}),dn(()=>{u.value=window.innerWidth,window.addEventListener("resize",Y),window.addEventListener("scroll",hs,{passive:!0})}),Tn(()=>{window.removeEventListener("resize",Y),window.removeEventListener("scroll",hs),d&&clearTimeout(d)}),(es,ms)=>{const _s=de("router-link");return E(),L("div",$y,[v("div",{class:"reading-progress",style:oa({width:r.value+"%"}),"aria-hidden":"true"},null,4),v("div",Iy,[v("div",Ny,[v("div",{class:"sticky-sidebar",ref_key:"leftSidebarContent",ref:m},[C.value?(E(),L("div",By,[bs(qr)])):rs("",!0)],512)]),v("div",zy,[v("div",qy,[v("div",Hy,[v("div",Vy,[w.value?(E(),L("div",Gy,[v("h1",Uy,G(w.value.title),1),v("div",Wy,[b.value&&w.value.date?(E(),L("span",Ky,[ms[0]||(ms[0]=v("i",{class:"fas fa-calendar-alt"},null,-1)),Ss(G(x.value)+" "+G(w.value.date),1)])):rs("",!0),Z.value>0?(E(),L("span",Yy,[ms[1]||(ms[1]=v("i",{class:"fas fa-clock"},null,-1)),Ss(G(B(Z.value)),1)])):rs("",!0),v("button",{type:"button",class:Ns(["article-copy-btn",{"article-copy-btn--copied":i.value}]),onClick:us,"aria-label":ls(n)("copyArticle"),"aria-live":"polite"},[v("i",{class:Ns(i.value?"fas fa-check":"fas fa-copy")},null,2),v("span",null,G(i.value?ls(n)("copied"):ls(n)("copyArticle")),1)],10,Xy)]),w.value.tags?.length?(E(),L("div",Jy,[(E(!0),L(js,null,Ts(w.value.tags,(fs,R)=>(E(),L("span",{key:R,class:"article-tag",onClick:H=>V(fs)}," # "+G(fs),9,Qy))),128))])):rs("",!0)])):rs("",!0),l.value?(E(),ma(yy,{key:1,rawMarkdown:l.value,articlePath:w.value?.path||"",articleTitle:w.value?.title||"",onMarkdownRendered:Vs},null,8,["rawMarkdown","articlePath","articleTitle"])):rs("",!0),l.value?(E(),L("nav",Zy,[!S.value&&F.value?(E(),L("span",sv)):rs("",!0),S.value?(E(),ma(_s,{key:1,to:O(S.value.path),class:"article-nav-item prev"},{default:hn(()=>[ms[2]||(ms[2]=v("div",{class:"nav-arrow"},[v("i",{class:"fas fa-arrow-left"})],-1)),v("div",nv,[v("div",av,G(f.value),1),v("div",ev,G(S.value.title),1)])]),_:1},8,["to"])):rs("",!0),F.value?(E(),ma(_s,{key:2,to:O(F.value.path),class:"article-nav-item next"},{default:hn(()=>[v("div",tv,[v("div",lv,G(k.value),1),v("div",ov,G(F.value.title),1)]),ms[3]||(ms[3]=v("div",{class:"nav-arrow"},[v("i",{class:"fas fa-arrow-right"})],-1))]),_:1},8,["to"])):rs("",!0)])):rs("",!0)])])])],512),v("div",pv,[v("div",{class:"sticky-sidebar",ref_key:"rightSidebarContent",ref:P},[C.value?(E(),L("div",cv,[bs(Vr,{ref_key:"onThisPageRef",ref:j,containerSelector:".markdown-body",levels:[2,3],offset:88},null,512)])):rs("",!0)],512)])]),l.value?(E(),ma(Fy,{key:0})):rs("",!0)])}}}),iv=mn(rv,[["__scopeId","data-v-bed67739"]]),Ba=()=>{if(typeof window<"u"){const s=localStorage.getItem("locale");if(s===_n[1])return Qn[1];if(s===_n[0])return Qn[0];if(navigator.language&&navigator.language.toLowerCase().startsWith("en"))return Qn[1]}return Qn[0]},Ep=s=>[{path:`${s}/`,name:`${s}-Home`,component:Ub},{path:`${s}/category`,name:`${s}-Category`,component:b_},{path:`${s}/resource`,name:`${s}-Resource`,component:q_},{path:`${s}/about`,name:`${s}-About`,component:fy},{path:`${s}/article/:path*`,name:`${s}-Article`,component:iv,props:!0}],uv=[{path:"/",redirect:()=>Ba()},{path:"/category",redirect:s=>`${Ba()}${s.path}`},{path:"/resource",redirect:s=>`${Ba()}${s.path}`},{path:"/about",redirect:s=>`${Ba()}${s.path}`},{path:"/article/:path*",redirect:s=>`${Ba()}${s.path}`},...Ep(Qn[0]),...Ep(Qn[1])],hv={mounted(s){if(typeof window>"u"||window.matchMedia("(prefers-reduced-motion: reduce)").matches||!("IntersectionObserver"in window))return;s.classList.add("reveal");const n=new IntersectionObserver(([a])=>{a.isIntersecting&&(s.classList.add("reveal-visible"),n.disconnect())},{threshold:.12});n.observe(s),s.__revealIo=n},unmounted(s){s.__revealIo?.disconnect(),s.__revealIo=null}};Mp(()=>import("./bootstrap.esm-inoyDHN5.js"),[]);am(Zf,{routes:uv},({app:s,router:n,initialState:a,isClient:e,routePath:t})=>{const l=em();s.use(l),a.pinia&&(l.state.value=a.pinia),s.use(n),s.use(Ye),s.mixin(Kh),s.directive("reveal",hv);const o=$l();o.initTheme(),e&&o.initLocale(),e&&"serviceWorker"in navigator&&window.addEventListener("load",()=>{navigator.serviceWorker.register("/sw.js").catch(()=>{})})});export{js as F,dv as T,Mp as _,fe as a,v as b,ma as c,zs as d,ls as e,Aa as f,L as g,Ss as h,rs as i,Ts as j,Lo as k,un as l,$n as m,Ks as n,E as o,mn as p,ns as r,G as t,Xs as u,mv as v,Bs as w};
