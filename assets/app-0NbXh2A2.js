const __vite__mapDeps=(i,m=__vite__mapDeps,d=(m.f||(m.f=["assets/SearchModal-BsN3_RSp.js","assets/SearchModal-CGQ99zei.css"])))=>i.map(i=>d[i]);
(function(){const n=document.createElement("link").relList;if(n&&n.supports&&n.supports("modulepreload"))return;for(const t of document.querySelectorAll('link[rel="modulepreload"]'))e(t);new MutationObserver(t=>{for(const l of t)if(l.type==="childList")for(const o of l.addedNodes)o.tagName==="LINK"&&o.rel==="modulepreload"&&e(o)}).observe(document,{childList:!0,subtree:!0});function a(t){const l={};return t.integrity&&(l.integrity=t.integrity),t.referrerPolicy&&(l.referrerPolicy=t.referrerPolicy),t.crossOrigin==="use-credentials"?l.credentials="include":t.crossOrigin==="anonymous"?l.credentials="omit":l.credentials="same-origin",l}function e(t){if(t.ep)return;t.ep=!0;const l=a(t);fetch(t.href,l)}})();const gu="modulepreload",fu=function(s){return"/"+s},Fo={},fs=function(n,a,e){let t=Promise.resolve();if(a&&a.length>0){let c=function(u){return Promise.all(u.map(r=>Promise.resolve(r).then(i=>({status:"fulfilled",value:i}),i=>({status:"rejected",reason:i}))))};document.getElementsByTagName("link");const o=document.querySelector("meta[property=csp-nonce]"),p=o?.nonce||o?.getAttribute("nonce");t=c(a.map(u=>{if(u=fu(u),u in Fo)return;Fo[u]=!0;const r=u.endsWith(".css"),i=r?'[rel="stylesheet"]':"";if(document.querySelector(`link[href="${u}"]${i}`))return;const h=document.createElement("link");if(h.rel=r?"stylesheet":gu,r||(h.as="script"),h.crossOrigin="",h.href=u,p&&h.setAttribute("nonce",p),document.head.appendChild(h),r)return new Promise((d,b)=>{h.addEventListener("load",d),h.addEventListener("error",()=>b(new Error(`Unable to preload CSS for ${u}`)))})}))}function l(o){const p=new Event("vite:preloadError",{cancelable:!0});if(p.payload=o,window.dispatchEvent(p),!p.defaultPrevented)throw o}return t.then(o=>{for(const p of o||[])p.status==="rejected"&&l(p.reason);return n().catch(l)})};/**
* @vue/shared v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/function Yl(s){const n=Object.create(null);for(const a of s.split(","))n[a]=1;return a=>a in n}const Ss={},Ua=[],$n=()=>{},Dc=()=>!1,St=s=>s.charCodeAt(0)===111&&s.charCodeAt(1)===110&&(s.charCodeAt(2)>122||s.charCodeAt(2)<97),Ql=s=>s.startsWith("onUpdate:"),Gs=Object.assign,Jl=(s,n)=>{const a=s.indexOf(n);a>-1&&s.splice(a,1)},bu=Object.prototype.hasOwnProperty,xs=(s,n)=>bu.call(s,n),is=Array.isArray,Ha=s=>Pt(s)==="[object Map]",Oc=s=>Pt(s)==="[object Set]",js=s=>typeof s=="function",Is=s=>typeof s=="string",ba=s=>typeof s=="symbol",Rs=s=>s!==null&&typeof s=="object",Ic=s=>(Rs(s)||js(s))&&js(s.then)&&js(s.catch),Nc=Object.prototype.toString,Pt=s=>Nc.call(s),_u=s=>Pt(s).slice(8,-1),Fc=s=>Pt(s)==="[object Object]",Zl=s=>Is(s)&&s!=="NaN"&&s[0]!=="-"&&""+parseInt(s,10)===s,ce=Yl(",key,ref,ref_for,ref_key,onVnodeBeforeMount,onVnodeMounted,onVnodeBeforeUpdate,onVnodeUpdated,onVnodeBeforeUnmount,onVnodeUnmounted"),Et=s=>{const n=Object.create(null);return(a=>n[a]||(n[a]=s(a)))},yu=/-\w/g,Pn=Et(s=>s.replace(yu,n=>n.slice(1).toUpperCase())),vu=/\B([A-Z])/g,_a=Et(s=>s.replace(vu,"-$1").toLowerCase()),Tt=Et(s=>s.charAt(0).toUpperCase()+s.slice(1)),Wt=Et(s=>s?`on${Tt(s)}`:""),da=(s,n)=>!Object.is(s,n),st=(s,...n)=>{for(let a=0;a<s.length;a++)s[a](...n)},$c=(s,n,a,e=!1)=>{Object.defineProperty(s,n,{configurable:!0,enumerable:!1,writable:e,value:a})},so=s=>{const n=parseFloat(s);return isNaN(n)?s:n},wu=s=>{const n=Is(s)?Number(s):NaN;return isNaN(n)?s:n};let $o;const At=()=>$o||($o=typeof globalThis<"u"?globalThis:typeof self<"u"?self:typeof window<"u"?window:typeof global<"u"?global:{});function qn(s){if(is(s)){const n={};for(let a=0;a<s.length;a++){const e=s[a],t=Is(e)?Su(e):qn(e);if(t)for(const l in t)n[l]=t[l]}return n}else if(Is(s)||Rs(s))return s}const Cu=/;(?![^(]*\))/g,ku=/:([^]+)/,xu=/\/\*[^]*?\*\//g;function Su(s){const n={};return s.replace(xu,"").split(Cu).forEach(a=>{if(a){const e=a.split(ku);e.length>1&&(n[e[0].trim()]=e[1].trim())}}),n}function Ts(s){let n="";if(Is(s))n=s;else if(is(s))for(let a=0;a<s.length;a++){const e=Ts(s[a]);e&&(n+=e+" ")}else if(Rs(s))for(const a in s)s[a]&&(n+=a+" ");return n.trim()}const Pu="itemscope,allowfullscreen,formnovalidate,ismap,nomodule,novalidate,readonly",Eu=Yl(Pu);function Bc(s){return!!s||s===""}const qc=s=>!!(s&&s.__v_isRef===!0),G=s=>Is(s)?s:s==null?"":is(s)||Rs(s)&&(s.toString===Nc||!js(s.toString))?qc(s)?G(s.value):JSON.stringify(s,zc,2):String(s),zc=(s,n)=>qc(n)?zc(s,n.value):Ha(n)?{[`Map(${n.size})`]:[...n.entries()].reduce((a,[e,t],l)=>(a[Kt(e,l)+" =>"]=t,a),{})}:Oc(n)?{[`Set(${n.size})`]:[...n.values()].map(a=>Kt(a))}:ba(n)?Kt(n):Rs(n)&&!is(n)&&!Fc(n)?String(n):n,Kt=(s,n="")=>{var a;return ba(s)?`Symbol(${(a=s.description)!=null?a:n})`:s};/**
* @vue/reactivity v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/let Zs;class Vc{constructor(n=!1){this.detached=n,this._active=!0,this._on=0,this.effects=[],this.cleanups=[],this._isPaused=!1,this.parent=Zs,!n&&Zs&&(this.index=(Zs.scopes||(Zs.scopes=[])).push(this)-1)}get active(){return this._active}pause(){if(this._active){this._isPaused=!0;let n,a;if(this.scopes)for(n=0,a=this.scopes.length;n<a;n++)this.scopes[n].pause();for(n=0,a=this.effects.length;n<a;n++)this.effects[n].pause()}}resume(){if(this._active&&this._isPaused){this._isPaused=!1;let n,a;if(this.scopes)for(n=0,a=this.scopes.length;n<a;n++)this.scopes[n].resume();for(n=0,a=this.effects.length;n<a;n++)this.effects[n].resume()}}run(n){if(this._active){const a=Zs;try{return Zs=this,n()}finally{Zs=a}}}on(){++this._on===1&&(this.prevScope=Zs,Zs=this)}off(){this._on>0&&--this._on===0&&(Zs=this.prevScope,this.prevScope=void 0)}stop(n){if(this._active){this._active=!1;let a,e;for(a=0,e=this.effects.length;a<e;a++)this.effects[a].stop();for(this.effects.length=0,a=0,e=this.cleanups.length;a<e;a++)this.cleanups[a]();if(this.cleanups.length=0,this.scopes){for(a=0,e=this.scopes.length;a<e;a++)this.scopes[a].stop(!0);this.scopes.length=0}if(!this.detached&&this.parent&&!n){const t=this.parent.scopes.pop();t&&t!==this&&(this.parent.scopes[this.index]=t,t.index=this.index)}this.parent=void 0}}}function no(s){return new Vc(s)}function Uc(){return Zs}function Tu(s,n=!1){Zs&&Zs.cleanups.push(s)}let Es;const Xt=new WeakSet;class Hc{constructor(n){this.fn=n,this.deps=void 0,this.depsTail=void 0,this.flags=5,this.next=void 0,this.cleanup=void 0,this.scheduler=void 0,Zs&&Zs.active&&Zs.effects.push(this)}pause(){this.flags|=64}resume(){this.flags&64&&(this.flags&=-65,Xt.has(this)&&(Xt.delete(this),this.trigger()))}notify(){this.flags&2&&!(this.flags&32)||this.flags&8||Wc(this)}run(){if(!(this.flags&1))return this.fn();this.flags|=2,Bo(this),Kc(this);const n=Es,a=An;Es=this,An=!0;try{return this.fn()}finally{Xc(this),Es=n,An=a,this.flags&=-3}}stop(){if(this.flags&1){for(let n=this.deps;n;n=n.nextDep)to(n);this.deps=this.depsTail=void 0,Bo(this),this.onStop&&this.onStop(),this.flags&=-2}}trigger(){this.flags&64?Xt.add(this):this.scheduler?this.scheduler():this.runIfDirty()}runIfDirty(){jl(this)&&this.run()}get dirty(){return jl(this)}}let Gc=0,re,ie;function Wc(s,n=!1){if(s.flags|=8,n){s.next=ie,ie=s;return}s.next=re,re=s}function ao(){Gc++}function eo(){if(--Gc>0)return;if(ie){let n=ie;for(ie=void 0;n;){const a=n.next;n.next=void 0,n.flags&=-9,n=a}}let s;for(;re;){let n=re;for(re=void 0;n;){const a=n.next;if(n.next=void 0,n.flags&=-9,n.flags&1)try{n.trigger()}catch(e){s||(s=e)}n=a}}if(s)throw s}function Kc(s){for(let n=s.deps;n;n=n.nextDep)n.version=-1,n.prevActiveLink=n.dep.activeLink,n.dep.activeLink=n}function Xc(s){let n,a=s.depsTail,e=a;for(;e;){const t=e.prevDep;e.version===-1?(e===a&&(a=t),to(e),Au(e)):n=e,e.dep.activeLink=e.prevActiveLink,e.prevActiveLink=void 0,e=t}s.deps=n,s.depsTail=a}function jl(s){for(let n=s.deps;n;n=n.nextDep)if(n.dep.version!==n.version||n.dep.computed&&(Yc(n.dep.computed)||n.dep.version!==n.version))return!0;return!!s._dirty}function Yc(s){if(s.flags&4&&!(s.flags&16)||(s.flags&=-17,s.globalVersion===ye)||(s.globalVersion=ye,!s.isSSR&&s.flags&128&&(!s.deps&&!s._dirty||!jl(s))))return;s.flags|=2;const n=s.dep,a=Es,e=An;Es=s,An=!0;try{Kc(s);const t=s.fn(s._value);(n.version===0||da(t,s._value))&&(s.flags|=128,s._value=t,n.version++)}catch(t){throw n.version++,t}finally{Es=a,An=e,Xc(s),s.flags&=-3}}function to(s,n=!1){const{dep:a,prevSub:e,nextSub:t}=s;if(e&&(e.nextSub=t,s.prevSub=void 0),t&&(t.prevSub=e,s.nextSub=void 0),a.subs===s&&(a.subs=e,!e&&a.computed)){a.computed.flags&=-5;for(let l=a.computed.deps;l;l=l.nextDep)to(l,!0)}!n&&!--a.sc&&a.map&&a.map.delete(a.key)}function Au(s){const{prevDep:n,nextDep:a}=s;n&&(n.nextDep=a,s.prevDep=void 0),a&&(a.prevDep=n,s.nextDep=void 0)}let An=!0;const Qc=[];function sa(){Qc.push(An),An=!1}function na(){const s=Qc.pop();An=s===void 0?!0:s}function Bo(s){const{cleanup:n}=s;if(s.cleanup=void 0,n){const a=Es;Es=void 0;try{n()}finally{Es=a}}}let ye=0;class Ru{constructor(n,a){this.sub=n,this.dep=a,this.version=a.version,this.nextDep=this.prevDep=this.nextSub=this.prevSub=this.prevActiveLink=void 0}}class lo{constructor(n){this.computed=n,this.version=0,this.activeLink=void 0,this.subs=void 0,this.map=void 0,this.key=void 0,this.sc=0,this.__v_skip=!0}track(n){if(!Es||!An||Es===this.computed)return;let a=this.activeLink;if(a===void 0||a.sub!==Es)a=this.activeLink=new Ru(Es,this),Es.deps?(a.prevDep=Es.depsTail,Es.depsTail.nextDep=a,Es.depsTail=a):Es.deps=Es.depsTail=a,Jc(a);else if(a.version===-1&&(a.version=this.version,a.nextDep)){const e=a.nextDep;e.prevDep=a.prevDep,a.prevDep&&(a.prevDep.nextDep=e),a.prevDep=Es.depsTail,a.nextDep=void 0,Es.depsTail.nextDep=a,Es.depsTail=a,Es.deps===a&&(Es.deps=e)}return a}trigger(n){this.version++,ye++,this.notify(n)}notify(n){ao();try{for(let a=this.subs;a;a=a.prevSub)a.sub.notify()&&a.sub.dep.notify()}finally{eo()}}}function Jc(s){if(s.dep.sc++,s.sub.flags&4){const n=s.dep.computed;if(n&&!s.dep.subs){n.flags|=20;for(let e=n.deps;e;e=e.nextDep)Jc(e)}const a=s.dep.subs;a!==s&&(s.prevSub=a,a&&(a.nextSub=s)),s.dep.subs=s}}const rt=new WeakMap,La=Symbol(""),gl=Symbol(""),ve=Symbol("");function nn(s,n,a){if(An&&Es){let e=rt.get(s);e||rt.set(s,e=new Map);let t=e.get(a);t||(e.set(a,t=new lo),t.map=e,t.key=a),t.track()}}function Yn(s,n,a,e,t,l){const o=rt.get(s);if(!o){ye++;return}const p=c=>{c&&c.trigger()};if(ao(),n==="clear")o.forEach(p);else{const c=is(s),u=c&&Zl(a);if(c&&a==="length"){const r=Number(e);o.forEach((i,h)=>{(h==="length"||h===ve||!ba(h)&&h>=r)&&p(i)})}else switch((a!==void 0||o.has(void 0))&&p(o.get(a)),u&&p(o.get(ve)),n){case"add":c?u&&p(o.get("length")):(p(o.get(La)),Ha(s)&&p(o.get(gl)));break;case"delete":c||(p(o.get(La)),Ha(s)&&p(o.get(gl)));break;case"set":Ha(s)&&p(o.get(La));break}}eo()}function Lu(s,n){const a=rt.get(s);return a&&a.get(n)}function Ia(s){const n=ws(s);return n===s?n:(nn(n,"iterate",ve),xn(s)?n:n.map(Ks))}function Rt(s){return nn(s=ws(s),"iterate",ve),s}const Mu={__proto__:null,[Symbol.iterator](){return Yt(this,Symbol.iterator,Ks)},concat(...s){return Ia(this).concat(...s.map(n=>is(n)?Ia(n):n))},entries(){return Yt(this,"entries",s=>(s[1]=Ks(s[1]),s))},every(s,n){return Un(this,"every",s,n,void 0,arguments)},filter(s,n){return Un(this,"filter",s,n,a=>a.map(Ks),arguments)},find(s,n){return Un(this,"find",s,n,Ks,arguments)},findIndex(s,n){return Un(this,"findIndex",s,n,void 0,arguments)},findLast(s,n){return Un(this,"findLast",s,n,Ks,arguments)},findLastIndex(s,n){return Un(this,"findLastIndex",s,n,void 0,arguments)},forEach(s,n){return Un(this,"forEach",s,n,void 0,arguments)},includes(...s){return Qt(this,"includes",s)},indexOf(...s){return Qt(this,"indexOf",s)},join(s){return Ia(this).join(s)},lastIndexOf(...s){return Qt(this,"lastIndexOf",s)},map(s,n){return Un(this,"map",s,n,void 0,arguments)},pop(){return ne(this,"pop")},push(...s){return ne(this,"push",s)},reduce(s,...n){return qo(this,"reduce",s,n)},reduceRight(s,...n){return qo(this,"reduceRight",s,n)},shift(){return ne(this,"shift")},some(s,n){return Un(this,"some",s,n,void 0,arguments)},splice(...s){return ne(this,"splice",s)},toReversed(){return Ia(this).toReversed()},toSorted(s){return Ia(this).toSorted(s)},toSpliced(...s){return Ia(this).toSpliced(...s)},unshift(...s){return ne(this,"unshift",s)},values(){return Yt(this,"values",Ks)}};function Yt(s,n,a){const e=Rt(s),t=e[n]();return e!==s&&!xn(s)&&(t._next=t.next,t.next=()=>{const l=t._next();return l.done||(l.value=a(l.value)),l}),t}const Du=Array.prototype;function Un(s,n,a,e,t,l){const o=Rt(s),p=o!==s&&!xn(s),c=o[n];if(c!==Du[n]){const i=c.apply(s,l);return p?Ks(i):i}let u=a;o!==s&&(p?u=function(i,h){return a.call(this,Ks(i),h,s)}:a.length>2&&(u=function(i,h){return a.call(this,i,h,s)}));const r=c.call(o,u,e);return p&&t?t(r):r}function qo(s,n,a,e){const t=Rt(s);let l=a;return t!==s&&(xn(s)?a.length>3&&(l=function(o,p,c){return a.call(this,o,p,c,s)}):l=function(o,p,c){return a.call(this,o,Ks(p),c,s)}),t[n](l,...e)}function Qt(s,n,a){const e=ws(s);nn(e,"iterate",ve);const t=e[n](...a);return(t===-1||t===!1)&&co(a[0])?(a[0]=ws(a[0]),e[n](...a)):t}function ne(s,n,a=[]){sa(),ao();const e=ws(s)[n].apply(s,a);return eo(),na(),e}const Ou=Yl("__proto__,__v_isRef,__isVue"),Zc=new Set(Object.getOwnPropertyNames(Symbol).filter(s=>s!=="arguments"&&s!=="caller").map(s=>Symbol[s]).filter(ba));function Iu(s){ba(s)||(s=String(s));const n=ws(this);return nn(n,"has",s),n.hasOwnProperty(s)}class sr{constructor(n=!1,a=!1){this._isReadonly=n,this._isShallow=a}get(n,a,e){if(a==="__v_skip")return n.__v_skip;const t=this._isReadonly,l=this._isShallow;if(a==="__v_isReactive")return!t;if(a==="__v_isReadonly")return t;if(a==="__v_isShallow")return l;if(a==="__v_raw")return e===(t?l?Gu:tr:l?er:ar).get(n)||Object.getPrototypeOf(n)===Object.getPrototypeOf(e)?n:void 0;const o=is(n);if(!t){let c;if(o&&(c=Mu[a]))return c;if(a==="hasOwnProperty")return Iu}const p=Reflect.get(n,a,Os(n)?n:e);if((ba(a)?Zc.has(a):Ou(a))||(t||nn(n,"get",a),l))return p;if(Os(p)){const c=o&&Zl(a)?p:p.value;return t&&Rs(c)?bl(c):c}return Rs(p)?t?bl(p):Me(p):p}}class nr extends sr{constructor(n=!1){super(!1,n)}set(n,a,e,t){let l=n[a];if(!this._isShallow){const c=ja(l);if(!xn(e)&&!ja(e)&&(l=ws(l),e=ws(e)),!is(n)&&Os(l)&&!Os(e))return c||(l.value=e),!0}const o=is(n)&&Zl(a)?Number(a)<n.length:xs(n,a),p=Reflect.set(n,a,e,Os(n)?n:t);return n===ws(t)&&(o?da(e,l)&&Yn(n,"set",a,e):Yn(n,"add",a,e)),p}deleteProperty(n,a){const e=xs(n,a);n[a];const t=Reflect.deleteProperty(n,a);return t&&e&&Yn(n,"delete",a,void 0),t}has(n,a){const e=Reflect.has(n,a);return(!ba(a)||!Zc.has(a))&&nn(n,"has",a),e}ownKeys(n){return nn(n,"iterate",is(n)?"length":La),Reflect.ownKeys(n)}}class Nu extends sr{constructor(n=!1){super(!0,n)}set(n,a){return!0}deleteProperty(n,a){return!0}}const Fu=new nr,$u=new Nu,Bu=new nr(!0);const fl=s=>s,Ue=s=>Reflect.getPrototypeOf(s);function qu(s,n,a){return function(...e){const t=this.__v_raw,l=ws(t),o=Ha(l),p=s==="entries"||s===Symbol.iterator&&o,c=s==="keys"&&o,u=t[s](...e),r=a?fl:n?it:Ks;return!n&&nn(l,"iterate",c?gl:La),{next(){const{value:i,done:h}=u.next();return h?{value:i,done:h}:{value:p?[r(i[0]),r(i[1])]:r(i),done:h}},[Symbol.iterator](){return this}}}}function He(s){return function(...n){return s==="delete"?!1:s==="clear"?void 0:this}}function zu(s,n){const a={get(t){const l=this.__v_raw,o=ws(l),p=ws(t);s||(da(t,p)&&nn(o,"get",t),nn(o,"get",p));const{has:c}=Ue(o),u=n?fl:s?it:Ks;if(c.call(o,t))return u(l.get(t));if(c.call(o,p))return u(l.get(p));l!==o&&l.get(t)},get size(){const t=this.__v_raw;return!s&&nn(ws(t),"iterate",La),t.size},has(t){const l=this.__v_raw,o=ws(l),p=ws(t);return s||(da(t,p)&&nn(o,"has",t),nn(o,"has",p)),t===p?l.has(t):l.has(t)||l.has(p)},forEach(t,l){const o=this,p=o.__v_raw,c=ws(p),u=n?fl:s?it:Ks;return!s&&nn(c,"iterate",La),p.forEach((r,i)=>t.call(l,u(r),u(i),o))}};return Gs(a,s?{add:He("add"),set:He("set"),delete:He("delete"),clear:He("clear")}:{add(t){!n&&!xn(t)&&!ja(t)&&(t=ws(t));const l=ws(this);return Ue(l).has.call(l,t)||(l.add(t),Yn(l,"add",t,t)),this},set(t,l){!n&&!xn(l)&&!ja(l)&&(l=ws(l));const o=ws(this),{has:p,get:c}=Ue(o);let u=p.call(o,t);u||(t=ws(t),u=p.call(o,t));const r=c.call(o,t);return o.set(t,l),u?da(l,r)&&Yn(o,"set",t,l):Yn(o,"add",t,l),this},delete(t){const l=ws(this),{has:o,get:p}=Ue(l);let c=o.call(l,t);c||(t=ws(t),c=o.call(l,t)),p&&p.call(l,t);const u=l.delete(t);return c&&Yn(l,"delete",t,void 0),u},clear(){const t=ws(this),l=t.size!==0,o=t.clear();return l&&Yn(t,"clear",void 0,void 0),o}}),["keys","values","entries",Symbol.iterator].forEach(t=>{a[t]=qu(t,s,n)}),a}function oo(s,n){const a=zu(s,n);return(e,t,l)=>t==="__v_isReactive"?!s:t==="__v_isReadonly"?s:t==="__v_raw"?e:Reflect.get(xs(a,t)&&t in e?a:e,t,l)}const Vu={get:oo(!1,!1)},Uu={get:oo(!1,!0)},Hu={get:oo(!0,!1)};const ar=new WeakMap,er=new WeakMap,tr=new WeakMap,Gu=new WeakMap;function Wu(s){switch(s){case"Object":case"Array":return 1;case"Map":case"Set":case"WeakMap":case"WeakSet":return 2;default:return 0}}function Ku(s){return s.__v_skip||!Object.isExtensible(s)?0:Wu(_u(s))}function Me(s){return ja(s)?s:po(s,!1,Fu,Vu,ar)}function lr(s){return po(s,!1,Bu,Uu,er)}function bl(s){return po(s,!0,$u,Hu,tr)}function po(s,n,a,e,t){if(!Rs(s)||s.__v_raw&&!(n&&s.__v_isReactive))return s;const l=Ku(s);if(l===0)return s;const o=t.get(s);if(o)return o;const p=new Proxy(s,l===2?e:a);return t.set(s,p),p}function ma(s){return ja(s)?ma(s.__v_raw):!!(s&&s.__v_isReactive)}function ja(s){return!!(s&&s.__v_isReadonly)}function xn(s){return!!(s&&s.__v_isShallow)}function co(s){return s?!!s.__v_raw:!1}function ws(s){const n=s&&s.__v_raw;return n?ws(n):s}function ro(s){return!xs(s,"__v_skip")&&Object.isExtensible(s)&&$c(s,"__v_skip",!0),s}const Ks=s=>Rs(s)?Me(s):s,it=s=>Rs(s)?bl(s):s;function Os(s){return s?s.__v_isRef===!0:!1}function rs(s){return or(s,!1)}function io(s){return or(s,!0)}function or(s,n){return Os(s)?s:new Xu(s,n)}class Xu{constructor(n,a){this.dep=new lo,this.__v_isRef=!0,this.__v_isShallow=!1,this._rawValue=a?n:ws(n),this._value=a?n:Ks(n),this.__v_isShallow=a}get value(){return this.dep.track(),this._value}set value(n){const a=this._rawValue,e=this.__v_isShallow||xn(n)||ja(n);n=e?n:ws(n),da(n,a)&&(this._rawValue=n,this._value=e?n:Ks(n),this.dep.trigger())}}function ss(s){return Os(s)?s.value:s}function Yu(s){return js(s)?s():ss(s)}const Qu={get:(s,n,a)=>n==="__v_raw"?s:ss(Reflect.get(s,n,a)),set:(s,n,a,e)=>{const t=s[n];return Os(t)&&!Os(a)?(t.value=a,!0):Reflect.set(s,n,a,e)}};function pr(s){return ma(s)?s:new Proxy(s,Qu)}function Ju(s){const n=is(s)?new Array(s.length):{};for(const a in s)n[a]=sh(s,a);return n}class Zu{constructor(n,a,e){this._object=n,this._key=a,this._defaultValue=e,this.__v_isRef=!0,this._value=void 0}get value(){const n=this._object[this._key];return this._value=n===void 0?this._defaultValue:n}set value(n){this._object[this._key]=n}get dep(){return Lu(ws(this._object),this._key)}}function sh(s,n,a){const e=s[n];return Os(e)?e:new Zu(s,n,a)}class nh{constructor(n,a,e){this.fn=n,this.setter=a,this._value=void 0,this.dep=new lo(this),this.__v_isRef=!0,this.deps=void 0,this.depsTail=void 0,this.flags=16,this.globalVersion=ye-1,this.next=void 0,this.effect=this,this.__v_isReadonly=!a,this.isSSR=e}notify(){if(this.flags|=16,!(this.flags&8)&&Es!==this)return Wc(this,!0),!0}get value(){const n=this.dep.track();return Yc(this),n&&(n.version=this.dep.version),this._value}set value(n){this.setter&&this.setter(n)}}function ah(s,n,a=!1){let e,t;return js(s)?e=s:(e=s.get,t=s.set),new nh(e,t,a)}const Ge={},ut=new WeakMap;let Ea;function eh(s,n=!1,a=Ea){if(a){let e=ut.get(a);e||ut.set(a,e=[]),e.push(s)}}function th(s,n,a=Ss){const{immediate:e,deep:t,once:l,scheduler:o,augmentJob:p,call:c}=a,u=y=>t?y:xn(y)||t===!1||t===0?Qn(y,1):Qn(y);let r,i,h,d,b=!1,f=!1;if(Os(s)?(i=()=>s.value,b=xn(s)):ma(s)?(i=()=>u(s),b=!0):is(s)?(f=!0,b=s.some(y=>ma(y)||xn(y)),i=()=>s.map(y=>{if(Os(y))return y.value;if(ma(y))return u(y);if(js(y))return c?c(y,2):y()})):js(s)?n?i=c?()=>c(s,2):s:i=()=>{if(h){sa();try{h()}finally{na()}}const y=Ea;Ea=r;try{return c?c(s,3,[d]):s(d)}finally{Ea=y}}:i=$n,n&&t){const y=i,x=t===!0?1/0:t;i=()=>Qn(y(),x)}const C=Uc(),m=()=>{r.stop(),C&&C.active&&Jl(C.effects,r)};if(l&&n){const y=n;n=(...x)=>{y(...x),m()}}let g=f?new Array(s.length).fill(Ge):Ge;const j=y=>{if(!(!(r.flags&1)||!r.dirty&&!y))if(n){const x=r.run();if(t||b||(f?x.some((R,E)=>da(R,g[E])):da(x,g))){h&&h();const R=Ea;Ea=r;try{const E=[x,g===Ge?void 0:f&&g[0]===Ge?[]:g,d];g=x,c?c(n,3,E):n(...E)}finally{Ea=R}}}else r.run()};return p&&p(j),r=new Hc(i),r.scheduler=o?()=>o(j,!1):j,d=y=>eh(y,!1,r),h=r.onStop=()=>{const y=ut.get(r);if(y){if(c)c(y,4);else for(const x of y)x();ut.delete(r)}},n?e?j(!0):g=r.run():o?o(j.bind(null,!0),!0):r.run(),m.pause=r.pause.bind(r),m.resume=r.resume.bind(r),m.stop=m,m}function Qn(s,n=1/0,a){if(n<=0||!Rs(s)||s.__v_skip||(a=a||new Map,(a.get(s)||0)>=n))return s;if(a.set(s,n),n--,Os(s))Qn(s.value,n,a);else if(is(s))for(let e=0;e<s.length;e++)Qn(s[e],n,a);else if(Oc(s)||Ha(s))s.forEach(e=>{Qn(e,n,a)});else if(Fc(s)){for(const e in s)Qn(s[e],n,a);for(const e of Object.getOwnPropertySymbols(s))Object.prototype.propertyIsEnumerable.call(s,e)&&Qn(s[e],n,a)}return s}/**
* @vue/runtime-core v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/function De(s,n,a,e){try{return e?s(...e):s()}catch(t){Oe(t,n,a)}}function Rn(s,n,a,e){if(js(s)){const t=De(s,n,a,e);return t&&Ic(t)&&t.catch(l=>{Oe(l,n,a)}),t}if(is(s)){const t=[];for(let l=0;l<s.length;l++)t.push(Rn(s[l],n,a,e));return t}}function Oe(s,n,a,e=!0){const t=n?n.vnode:null,{errorHandler:l,throwUnhandledErrorInProduction:o}=n&&n.appContext.config||Ss;if(n){let p=n.parent;const c=n.proxy,u=`https://vuejs.org/error-reference/#runtime-${a}`;for(;p;){const r=p.ec;if(r){for(let i=0;i<r.length;i++)if(r[i](s,c,u)===!1)return}p=p.parent}if(l){sa(),De(l,null,10,[s,c,u]),na();return}}lh(s,a,t,e,o)}function lh(s,n,a,e=!0,t=!1){if(t)throw s;console.error(s)}const rn=[];let Nn=-1;const Ga=[];let ia=null,Ba=0;const cr=Promise.resolve();let ht=null;function bn(s){const n=ht||cr;return s?n.then(this?s.bind(this):s):n}function oh(s){let n=Nn+1,a=rn.length;for(;n<a;){const e=n+a>>>1,t=rn[e],l=we(t);l<s||l===s&&t.flags&2?n=e+1:a=e}return n}function uo(s){if(!(s.flags&1)){const n=we(s),a=rn[rn.length-1];!a||!(s.flags&2)&&n>=we(a)?rn.push(s):rn.splice(oh(n),0,s),s.flags|=1,rr()}}function rr(){ht||(ht=cr.then(ur))}function ph(s){is(s)?Ga.push(...s):ia&&s.id===-1?ia.splice(Ba+1,0,s):s.flags&1||(Ga.push(s),s.flags|=1),rr()}function zo(s,n,a=Nn+1){for(;a<rn.length;a++){const e=rn[a];if(e&&e.flags&2){if(s&&e.id!==s.uid)continue;rn.splice(a,1),a--,e.flags&4&&(e.flags&=-2),e(),e.flags&4||(e.flags&=-2)}}}function ir(s){if(Ga.length){const n=[...new Set(Ga)].sort((a,e)=>we(a)-we(e));if(Ga.length=0,ia){ia.push(...n);return}for(ia=n,Ba=0;Ba<ia.length;Ba++){const a=ia[Ba];a.flags&4&&(a.flags&=-2),a.flags&8||a(),a.flags&=-2}ia=null,Ba=0}}const we=s=>s.id==null?s.flags&2?-1:1/0:s.id;function ur(s){try{for(Nn=0;Nn<rn.length;Nn++){const n=rn[Nn];n&&!(n.flags&8)&&(n.flags&4&&(n.flags&=-2),De(n,n.i,n.i?15:14),n.flags&4||(n.flags&=-2))}}finally{for(;Nn<rn.length;Nn++){const n=rn[Nn];n&&(n.flags&=-2)}Nn=-1,rn.length=0,ir(),ht=null,(rn.length||Ga.length)&&ur()}}let _n=null,hr=null;function dt(s){const n=_n;return _n=s,hr=s&&s.type.__scopeId||null,n}function jn(s,n=_n,a){if(!n||s._n)return s;const e=(...t)=>{e._d&&gt(-1);const l=dt(n);let o;try{o=s(...t)}finally{dt(l),e._d&&gt(1)}return o};return e._n=!0,e._c=!0,e._d=!0,e}function an(s,n){if(_n===null)return s;const a=Dt(_n),e=s.dirs||(s.dirs=[]);for(let t=0;t<n.length;t++){let[l,o,p,c=Ss]=n[t];l&&(js(l)&&(l={mounted:l,updated:l}),l.deep&&Qn(o),e.push({dir:l,instance:a,value:o,oldValue:void 0,arg:p,modifiers:c}))}return s}function ka(s,n,a,e){const t=s.dirs,l=n&&n.dirs;for(let o=0;o<t.length;o++){const p=t[o];l&&(p.oldValue=l[o].value);let c=p.dir[e];c&&(sa(),Rn(c,a,8,[s.el,p,s,n]),na())}}const dr=Symbol("_vte"),mr=s=>s.__isTeleport,ue=s=>s&&(s.disabled||s.disabled===""),Vo=s=>s&&(s.defer||s.defer===""),Uo=s=>typeof SVGElement<"u"&&s instanceof SVGElement,Ho=s=>typeof MathMLElement=="function"&&s instanceof MathMLElement,_l=(s,n)=>{const a=s&&s.to;return Is(a)?n?n(a):null:a},jr={name:"Teleport",__isTeleport:!0,process(s,n,a,e,t,l,o,p,c,u){const{mc:r,pc:i,pbc:h,o:{insert:d,querySelector:b,createText:f,createComment:C}}=u,m=ue(n.props);let{shapeFlag:g,children:j,dynamicChildren:y}=n;if(s==null){const x=n.el=f(""),R=n.anchor=f("");d(x,a,e),d(R,a,e);const E=(P,U)=>{g&16&&r(j,P,U,t,l,o,p,c)},O=()=>{const P=n.target=_l(n.props,b),U=gr(P,n,f,d);P&&(o!=="svg"&&Uo(P)?o="svg":o!=="mathml"&&Ho(P)&&(o="mathml"),t&&t.isCE&&(t.ce._teleportTargets||(t.ce._teleportTargets=new Set)).add(P),m||(E(P,U),nt(n,!1)))};m&&(E(a,R),nt(n,!0)),Vo(n.props)?(n.el.__isMounted=!1,pn(()=>{O(),delete n.el.__isMounted},l)):O()}else{if(Vo(n.props)&&s.el.__isMounted===!1){pn(()=>{jr.process(s,n,a,e,t,l,o,p,c,u)},l);return}n.el=s.el,n.targetStart=s.targetStart;const x=n.anchor=s.anchor,R=n.target=s.target,E=n.targetAnchor=s.targetAnchor,O=ue(s.props),P=O?a:R,U=O?x:E;if(o==="svg"||Uo(R)?o="svg":(o==="mathml"||Ho(R))&&(o="mathml"),y?(h(s.dynamicChildren,y,P,t,l,o,p),_o(s,n,!0)):c||i(s,n,P,U,t,l,o,p,!1),m)O?n.props&&s.props&&n.props.to!==s.props.to&&(n.props.to=s.props.to):We(n,a,x,u,1);else if((n.props&&n.props.to)!==(s.props&&s.props.to)){const Y=n.target=_l(n.props,b);Y&&We(n,Y,null,u,0)}else O&&We(n,R,E,u,1);nt(n,m)}},remove(s,n,a,{um:e,o:{remove:t}},l){const{shapeFlag:o,children:p,anchor:c,targetStart:u,targetAnchor:r,target:i,props:h}=s;if(i&&(t(u),t(r)),l&&t(c),o&16){const d=l||!ue(h);for(let b=0;b<p.length;b++){const f=p[b];e(f,n,a,d,!!f.dynamicChildren)}}},move:We,hydrate:ch};function We(s,n,a,{o:{insert:e},m:t},l=2){l===0&&e(s.targetAnchor,n,a);const{el:o,anchor:p,shapeFlag:c,children:u,props:r}=s,i=l===2;if(i&&e(o,n,a),(!i||ue(r))&&c&16)for(let h=0;h<u.length;h++)t(u[h],n,a,2);i&&e(p,n,a)}function ch(s,n,a,e,t,l,{o:{nextSibling:o,parentNode:p,querySelector:c,insert:u,createText:r}},i){function h(f,C,m,g){C.anchor=i(o(f),C,p(f),a,e,t,l),C.targetStart=m,C.targetAnchor=g}const d=n.target=_l(n.props,c),b=ue(n.props);if(d){const f=d._lpa||d.firstChild;if(n.shapeFlag&16)if(b)h(s,n,f,f&&o(f));else{n.anchor=o(s);let C=f;for(;C;){if(C&&C.nodeType===8){if(C.data==="teleport start anchor")n.targetStart=C;else if(C.data==="teleport anchor"){n.targetAnchor=C,d._lpa=n.targetAnchor&&o(n.targetAnchor);break}}C=o(C)}n.targetAnchor||gr(d,n,r,u),i(f&&o(f),n,d,a,e,t,l)}nt(n,b)}else b&&n.shapeFlag&16&&h(s,n,s,o(s));return n.anchor&&o(n.anchor)}const V1=jr;function nt(s,n){const a=s.ctx;if(a&&a.ut){let e,t;for(n?(e=s.el,t=s.anchor):(e=s.targetStart,t=s.targetAnchor);e&&e!==t;)e.nodeType===1&&e.setAttribute("data-v-owner",a.uid),e=e.nextSibling;a.ut()}}function gr(s,n,a,e){const t=n.targetStart=a(""),l=n.targetAnchor=a("");return t[dr]=l,s&&(e(t,s),e(l,s)),l}const Xn=Symbol("_leaveCb"),Ke=Symbol("_enterCb");function rh(){const s={isMounted:!1,isLeaving:!1,isUnmounting:!1,leavingVNodes:new Map};return En(()=>{s.isMounted=!0}),zn(()=>{s.isUnmounting=!0}),s}const Cn=[Function,Array],fr={mode:String,appear:Boolean,persisted:Boolean,onBeforeEnter:Cn,onEnter:Cn,onAfterEnter:Cn,onEnterCancelled:Cn,onBeforeLeave:Cn,onLeave:Cn,onAfterLeave:Cn,onLeaveCancelled:Cn,onBeforeAppear:Cn,onAppear:Cn,onAfterAppear:Cn,onAppearCancelled:Cn},br=s=>{const n=s.subTree;return n.component?br(n.component):n},ih={name:"BaseTransition",props:fr,setup(s,{slots:n}){const a=la(),e=rh();return()=>{const t=n.default&&vr(n.default(),!0);if(!t||!t.length)return;const l=_r(t),o=ws(s),{mode:p}=o;if(e.isLeaving)return Jt(l);const c=Go(l);if(!c)return Jt(l);let u=yl(c,o,e,a,i=>u=i);c.type!==un&&Ce(c,u);let r=a.subTree&&Go(a.subTree);if(r&&r.type!==un&&!Ta(r,c)&&br(a).type!==un){let i=yl(r,o,e,a);if(Ce(r,i),p==="out-in"&&c.type!==un)return e.isLeaving=!0,i.afterLeave=()=>{e.isLeaving=!1,a.job.flags&8||a.update(),delete i.afterLeave,r=void 0},Jt(l);p==="in-out"&&c.type!==un?i.delayLeave=(h,d,b)=>{const f=yr(e,r);f[String(r.key)]=r,h[Xn]=()=>{d(),h[Xn]=void 0,delete u.delayedLeave,r=void 0},u.delayedLeave=()=>{b(),delete u.delayedLeave,r=void 0}}:r=void 0}else r&&(r=void 0);return l}}};function _r(s){let n=s[0];if(s.length>1){for(const a of s)if(a.type!==un){n=a;break}}return n}const uh=ih;function yr(s,n){const{leavingVNodes:a}=s;let e=a.get(n.type);return e||(e=Object.create(null),a.set(n.type,e)),e}function yl(s,n,a,e,t){const{appear:l,mode:o,persisted:p=!1,onBeforeEnter:c,onEnter:u,onAfterEnter:r,onEnterCancelled:i,onBeforeLeave:h,onLeave:d,onAfterLeave:b,onLeaveCancelled:f,onBeforeAppear:C,onAppear:m,onAfterAppear:g,onAppearCancelled:j}=n,y=String(s.key),x=yr(a,s),R=(P,U)=>{P&&Rn(P,e,9,U)},E=(P,U)=>{const Y=U[1];R(P,U),is(P)?P.every(V=>V.length<=1)&&Y():P.length<=1&&Y()},O={mode:o,persisted:p,beforeEnter(P){let U=c;if(!a.isMounted)if(l)U=C||c;else return;P[Xn]&&P[Xn](!0);const Y=x[y];Y&&Ta(s,Y)&&Y.el[Xn]&&Y.el[Xn](),R(U,[P])},enter(P){let U=u,Y=r,V=i;if(!a.isMounted)if(l)U=m||u,Y=g||r,V=j||i;else return;let z=!1;const Q=P[Ke]=cs=>{z||(z=!0,cs?R(V,[P]):R(Y,[P]),O.delayedLeave&&O.delayedLeave(),P[Ke]=void 0)};U?E(U,[P,Q]):Q()},leave(P,U){const Y=String(s.key);if(P[Ke]&&P[Ke](!0),a.isUnmounting)return U();R(h,[P]);let V=!1;const z=P[Xn]=Q=>{V||(V=!0,U(),Q?R(f,[P]):R(b,[P]),P[Xn]=void 0,x[Y]===s&&delete x[Y])};x[Y]=s,d?E(d,[P,z]):z()},clone(P){const U=yl(P,n,a,e,t);return t&&t(U),U}};return O}function Jt(s){if(Ie(s))return s=ga(s),s.children=null,s}function Go(s){if(!Ie(s))return mr(s.type)&&s.children?_r(s.children):s;if(s.component)return s.component.subTree;const{shapeFlag:n,children:a}=s;if(a){if(n&16)return a[0];if(n&32&&js(a.default))return a.default()}}function Ce(s,n){s.shapeFlag&6&&s.component?(s.transition=n,Ce(s.component.subTree,n)):s.shapeFlag&128?(s.ssContent.transition=n.clone(s.ssContent),s.ssFallback.transition=n.clone(s.ssFallback)):s.transition=n}function vr(s,n=!1,a){let e=[],t=0;for(let l=0;l<s.length;l++){let o=s[l];const p=a==null?o.key:String(a)+String(o.key!=null?o.key:l);o.type===gs?(o.patchFlag&128&&t++,e=e.concat(vr(o.children,n,p))):(n||o.type!==un)&&e.push(p!=null?ga(o,{key:p}):o)}if(t>1)for(let l=0;l<e.length;l++)e[l].patchFlag=-2;return e}function Ns(s,n){return js(s)?Gs({name:s.name},n,{setup:s}):s}function ho(s){s.ids=[s.ids[0]+s.ids[2]+++"-",0,0]}function at(s){const n=la(),a=io(null);if(n){const t=n.refs===Ss?n.refs={}:n.refs;Object.defineProperty(t,s,{enumerable:!0,get:()=>a.value,set:l=>a.value=l})}return a}const mt=new WeakMap;function he(s,n,a,e,t=!1){if(is(s)){s.forEach((b,f)=>he(b,n&&(is(n)?n[f]:n),a,e,t));return}if(de(e)&&!t){e.shapeFlag&512&&e.type.__asyncResolved&&e.component.subTree.component&&he(s,n,a,e.component.subTree);return}const l=e.shapeFlag&4?Dt(e.component):e.el,o=t?null:l,{i:p,r:c}=s,u=n&&n.r,r=p.refs===Ss?p.refs={}:p.refs,i=p.setupState,h=ws(i),d=i===Ss?Dc:b=>xs(h,b);if(u!=null&&u!==c){if(Wo(n),Is(u))r[u]=null,d(u)&&(i[u]=null);else if(Os(u)){u.value=null;const b=n;b.k&&(r[b.k]=null)}}if(js(c))De(c,p,12,[o,r]);else{const b=Is(c),f=Os(c);if(b||f){const C=()=>{if(s.f){const m=b?d(c)?i[c]:r[c]:c.value;if(t)is(m)&&Jl(m,l);else if(is(m))m.includes(l)||m.push(l);else if(b)r[c]=[l],d(c)&&(i[c]=r[c]);else{const g=[l];c.value=g,s.k&&(r[s.k]=g)}}else b?(r[c]=o,d(c)&&(i[c]=o)):f&&(c.value=o,s.k&&(r[s.k]=o))};if(o){const m=()=>{C(),mt.delete(s)};m.id=-1,mt.set(s,m),pn(m,a)}else Wo(s),C()}}}function Wo(s){const n=mt.get(s);n&&(n.flags|=8,mt.delete(s))}const Ko=s=>s.nodeType===8;At().requestIdleCallback;At().cancelIdleCallback;function hh(s,n){if(Ko(s)&&s.data==="["){let a=1,e=s.nextSibling;for(;e;){if(e.nodeType===1){if(n(e)===!1)break}else if(Ko(e))if(e.data==="]"){if(--a===0)break}else e.data==="["&&a++;e=e.nextSibling}}else n(s)}const de=s=>!!s.type.__asyncLoader;function dh(s){js(s)&&(s={loader:s});const{loader:n,loadingComponent:a,errorComponent:e,delay:t=200,hydrate:l,timeout:o,suspensible:p=!0,onError:c}=s;let u=null,r,i=0;const h=()=>(i++,u=null,d()),d=()=>{let b;return u||(b=u=n().catch(f=>{if(f=f instanceof Error?f:new Error(String(f)),c)return new Promise((C,m)=>{c(f,()=>C(h()),()=>m(f),i+1)});throw f}).then(f=>b!==u&&u?u:(f&&(f.__esModule||f[Symbol.toStringTag]==="Module")&&(f=f.default),r=f,f)))};return Ns({name:"AsyncComponentWrapper",__asyncLoader:d,__asyncHydrate(b,f,C){let m=!1;(f.bu||(f.bu=[])).push(()=>m=!0);const g=()=>{m||C()},j=l?()=>{const y=l(g,x=>hh(b,x));y&&(f.bum||(f.bum=[])).push(y)}:g;r?j():d().then(()=>!f.isUnmounted&&j())},get __asyncResolved(){return r},setup(){const b=Xs;if(ho(b),r)return()=>Xe(r,b);const f=j=>{u=null,Oe(j,b,13,!e)};if(p&&b.suspense||Ka)return d().then(j=>()=>Xe(j,b)).catch(j=>(f(j),()=>e?bs(e,{error:j}):null));const C=rs(!1),m=rs(),g=rs(!!t);return t&&setTimeout(()=>{g.value=!1},t),o!=null&&setTimeout(()=>{if(!C.value&&!m.value){const j=new Error(`Async component timed out after ${o}ms.`);f(j),m.value=j}},o),d().then(()=>{C.value=!0,b.parent&&Ie(b.parent.vnode)&&b.parent.update()}).catch(j=>{f(j),m.value=j}),()=>{if(C.value&&r)return Xe(r,b);if(m.value&&e)return bs(e,{error:m.value});if(a&&!g.value)return Xe(a,b)}}})}function Xe(s,n){const{ref:a,props:e,children:t,ce:l}=n.vnode,o=bs(s,e,t);return o.ref=a,o.ce=l,delete n.vnode.ce,o}const Ie=s=>s.type.__isKeepAlive;function wr(s,n){kr(s,"a",n)}function Cr(s,n){kr(s,"da",n)}function kr(s,n,a=Xs){const e=s.__wdc||(s.__wdc=()=>{let t=a;for(;t;){if(t.isDeactivated)return;t=t.parent}return s()});if(Lt(n,e,a),a){let t=a.parent;for(;t&&t.parent;)Ie(t.parent.vnode)&&mh(e,n,a,t),t=t.parent}}function mh(s,n,a,e){const t=Lt(n,s,e,!0);mo(()=>{Jl(e[n],t)},a)}function Lt(s,n,a=Xs,e=!1){if(a){const t=a[s]||(a[s]=[]),l=n.__weh||(n.__weh=(...o)=>{sa();const p=$e(a),c=Rn(n,a,s,o);return p(),na(),c});return e?t.unshift(l):t.push(l),l}}const ta=s=>(n,a=Xs)=>{(!Ka||s==="sp")&&Lt(s,(...e)=>n(...e),a)},jh=ta("bm"),En=ta("m"),gh=ta("bu"),fh=ta("u"),zn=ta("bum"),mo=ta("um"),bh=ta("sp"),_h=ta("rtg"),yh=ta("rtc");function vh(s,n=Xs){Lt("ec",s,n)}const jo="components",wh="directives";function Za(s,n){return go(jo,s,!0,n)||s}const xr=Symbol.for("v-ndc");function Ch(s){return Is(s)?go(jo,s,!1)||s:s||xr}function Ne(s){return go(wh,s)}function go(s,n,a=!0,e=!1){const t=_n||Xs;if(t){const l=t.type;if(s===jo){const p=ud(l,!1);if(p&&(p===n||p===Pn(n)||p===Tt(Pn(n))))return l}const o=Xo(t[s]||l[s],n)||Xo(t.appContext[s],n);return!o&&e?l:o}}function Xo(s,n){return s&&(s[n]||s[Pn(n)]||s[Tt(Pn(n))])}function As(s,n,a,e){let t;const l=a,o=is(s);if(o||Is(s)){const p=o&&ma(s);let c=!1,u=!1;p&&(c=!xn(s),u=ja(s),s=Rt(s)),t=new Array(s.length);for(let r=0,i=s.length;r<i;r++)t[r]=n(c?u?it(Ks(s[r])):Ks(s[r]):s[r],r,void 0,l)}else if(typeof s=="number"){t=new Array(s);for(let p=0;p<s;p++)t[p]=n(p+1,p,void 0,l)}else if(Rs(s))if(s[Symbol.iterator])t=Array.from(s,(p,c)=>n(p,c,void 0,l));else{const p=Object.keys(s);t=new Array(p.length);for(let c=0,u=p.length;c<u;c++){const r=p[c];t[c]=n(s[r],r,c,l)}}else t=[];return t}const vl=s=>s?Ur(s)?Dt(s):vl(s.parent):null,me=Gs(Object.create(null),{$:s=>s,$el:s=>s.vnode.el,$data:s=>s.data,$props:s=>s.props,$attrs:s=>s.attrs,$slots:s=>s.slots,$refs:s=>s.refs,$parent:s=>vl(s.parent),$root:s=>vl(s.root),$host:s=>s.ce,$emit:s=>s.emit,$options:s=>Pr(s),$forceUpdate:s=>s.f||(s.f=()=>{uo(s.update)}),$nextTick:s=>s.n||(s.n=bn.bind(s.proxy)),$watch:s=>Hh.bind(s)}),Zt=(s,n)=>s!==Ss&&!s.__isScriptSetup&&xs(s,n),kh={get({_:s},n){if(n==="__v_skip")return!0;const{ctx:a,setupState:e,data:t,props:l,accessCache:o,type:p,appContext:c}=s;let u;if(n[0]!=="$"){const d=o[n];if(d!==void 0)switch(d){case 1:return e[n];case 2:return t[n];case 4:return a[n];case 3:return l[n]}else{if(Zt(e,n))return o[n]=1,e[n];if(t!==Ss&&xs(t,n))return o[n]=2,t[n];if((u=s.propsOptions[0])&&xs(u,n))return o[n]=3,l[n];if(a!==Ss&&xs(a,n))return o[n]=4,a[n];wl&&(o[n]=0)}}const r=me[n];let i,h;if(r)return n==="$attrs"&&nn(s.attrs,"get",""),r(s);if((i=p.__cssModules)&&(i=i[n]))return i;if(a!==Ss&&xs(a,n))return o[n]=4,a[n];if(h=c.config.globalProperties,xs(h,n))return h[n]},set({_:s},n,a){const{data:e,setupState:t,ctx:l}=s;return Zt(t,n)?(t[n]=a,!0):e!==Ss&&xs(e,n)?(e[n]=a,!0):xs(s.props,n)||n[0]==="$"&&n.slice(1)in s?!1:(l[n]=a,!0)},has({_:{data:s,setupState:n,accessCache:a,ctx:e,appContext:t,propsOptions:l,type:o}},p){let c,u;return!!(a[p]||s!==Ss&&p[0]!=="$"&&xs(s,p)||Zt(n,p)||(c=l[0])&&xs(c,p)||xs(e,p)||xs(me,p)||xs(t.config.globalProperties,p)||(u=o.__cssModules)&&u[p])},defineProperty(s,n,a){return a.get!=null?s._.accessCache[n]=0:xs(a,"value")&&this.set(s,n,a.value,null),Reflect.defineProperty(s,n,a)}};function Yo(s){return is(s)?s.reduce((n,a)=>(n[a]=null,n),{}):s}let wl=!0;function xh(s){const n=Pr(s),a=s.proxy,e=s.ctx;wl=!1,n.beforeCreate&&Qo(n.beforeCreate,s,"bc");const{data:t,computed:l,methods:o,watch:p,provide:c,inject:u,created:r,beforeMount:i,mounted:h,beforeUpdate:d,updated:b,activated:f,deactivated:C,beforeDestroy:m,beforeUnmount:g,destroyed:j,unmounted:y,render:x,renderTracked:R,renderTriggered:E,errorCaptured:O,serverPrefetch:P,expose:U,inheritAttrs:Y,components:V,directives:z,filters:Q}=n;if(u&&Sh(u,e,null),o)for(const X in o){const us=o[X];js(us)&&(e[X]=us.bind(a))}if(t){const X=t.call(a,a);Rs(X)&&(s.data=Me(X))}if(wl=!0,l)for(const X in l){const us=l[X],Fs=js(us)?us.bind(a,a):js(us.get)?us.get.bind(a,a):$n,Qs=!js(us)&&js(us.set)?us.set.bind(a):$n,Ms=ns({get:Fs,set:Qs});Object.defineProperty(e,X,{enumerable:!0,configurable:!0,get:()=>Ms.value,set:$s=>Ms.value=$s})}if(p)for(const X in p)Sr(p[X],e,a,X);if(c){const X=js(c)?c.call(a):c;Reflect.ownKeys(X).forEach(us=>{et(us,X[us])})}r&&Qo(r,s,"c");function ls(X,us){is(us)?us.forEach(Fs=>X(Fs.bind(a))):us&&X(us.bind(a))}if(ls(jh,i),ls(En,h),ls(gh,d),ls(fh,b),ls(wr,f),ls(Cr,C),ls(vh,O),ls(yh,R),ls(_h,E),ls(zn,g),ls(mo,y),ls(bh,P),is(U))if(U.length){const X=s.exposed||(s.exposed={});U.forEach(us=>{Object.defineProperty(X,us,{get:()=>a[us],set:Fs=>a[us]=Fs,enumerable:!0})})}else s.exposed||(s.exposed={});x&&s.render===$n&&(s.render=x),Y!=null&&(s.inheritAttrs=Y),V&&(s.components=V),z&&(s.directives=z),P&&ho(s)}function Sh(s,n,a=$n){is(s)&&(s=Cl(s));for(const e in s){const t=s[e];let l;Rs(t)?"default"in t?l=gn(t.from||e,t.default,!0):l=gn(t.from||e):l=gn(t),Os(l)?Object.defineProperty(n,e,{enumerable:!0,configurable:!0,get:()=>l.value,set:o=>l.value=o}):n[e]=l}}function Qo(s,n,a){Rn(is(s)?s.map(e=>e.bind(n.proxy)):s.bind(n.proxy),n,a)}function Sr(s,n,a,e){let t=e.includes(".")?$r(a,e):()=>a[e];if(Is(s)){const l=n[s];js(l)&&Vs(t,l)}else if(js(s))Vs(t,s.bind(a));else if(Rs(s))if(is(s))s.forEach(l=>Sr(l,n,a,e));else{const l=js(s.handler)?s.handler.bind(a):n[s.handler];js(l)&&Vs(t,l,s)}}function Pr(s){const n=s.type,{mixins:a,extends:e}=n,{mixins:t,optionsCache:l,config:{optionMergeStrategies:o}}=s.appContext,p=l.get(n);let c;return p?c=p:!t.length&&!a&&!e?c=n:(c={},t.length&&t.forEach(u=>jt(c,u,o,!0)),jt(c,n,o)),Rs(n)&&l.set(n,c),c}function jt(s,n,a,e=!1){const{mixins:t,extends:l}=n;l&&jt(s,l,a,!0),t&&t.forEach(o=>jt(s,o,a,!0));for(const o in n)if(!(e&&o==="expose")){const p=Ph[o]||a&&a[o];s[o]=p?p(s[o],n[o]):n[o]}return s}const Ph={data:Jo,props:Zo,emits:Zo,methods:pe,computed:pe,beforeCreate:ln,created:ln,beforeMount:ln,mounted:ln,beforeUpdate:ln,updated:ln,beforeDestroy:ln,beforeUnmount:ln,destroyed:ln,unmounted:ln,activated:ln,deactivated:ln,errorCaptured:ln,serverPrefetch:ln,components:pe,directives:pe,watch:Th,provide:Jo,inject:Eh};function Jo(s,n){return n?s?function(){return Gs(js(s)?s.call(this,this):s,js(n)?n.call(this,this):n)}:n:s}function Eh(s,n){return pe(Cl(s),Cl(n))}function Cl(s){if(is(s)){const n={};for(let a=0;a<s.length;a++)n[s[a]]=s[a];return n}return s}function ln(s,n){return s?[...new Set([].concat(s,n))]:n}function pe(s,n){return s?Gs(Object.create(null),s,n):n}function Zo(s,n){return s?is(s)&&is(n)?[...new Set([...s,...n])]:Gs(Object.create(null),Yo(s),Yo(n??{})):n}function Th(s,n){if(!s)return n;if(!n)return s;const a=Gs(Object.create(null),s);for(const e in n)a[e]=ln(s[e],n[e]);return a}function Er(){return{app:null,config:{isNativeTag:Dc,performance:!1,globalProperties:{},optionMergeStrategies:{},errorHandler:void 0,warnHandler:void 0,compilerOptions:{}},mixins:[],components:{},directives:{},provides:Object.create(null),optionsCache:new WeakMap,propsCache:new WeakMap,emitsCache:new WeakMap}}let Ah=0;function Rh(s,n){return function(e,t=null){js(e)||(e=Gs({},e)),t!=null&&!Rs(t)&&(t=null);const l=Er(),o=new WeakSet,p=[];let c=!1;const u=l.app={_uid:Ah++,_component:e,_props:t,_container:null,_context:l,_instance:null,version:dd,get config(){return l.config},set config(r){},use(r,...i){return o.has(r)||(r&&js(r.install)?(o.add(r),r.install(u,...i)):js(r)&&(o.add(r),r(u,...i))),u},mixin(r){return l.mixins.includes(r)||l.mixins.push(r),u},component(r,i){return i?(l.components[r]=i,u):l.components[r]},directive(r,i){return i?(l.directives[r]=i,u):l.directives[r]},mount(r,i,h){if(!c){const d=u._ceVNode||bs(e,t);return d.appContext=l,h===!0?h="svg":h===!1&&(h=void 0),s(d,r,h),c=!0,u._container=r,r.__vue_app__=u,Dt(d.component)}},onUnmount(r){p.push(r)},unmount(){c&&(Rn(p,u._instance,16),s(null,u._container),delete u._container.__vue_app__)},provide(r,i){return l.provides[r]=i,u},runWithContext(r){const i=Ma;Ma=u;try{return r()}finally{Ma=i}}};return u}}let Ma=null;function et(s,n){if(Xs){let a=Xs.provides;const e=Xs.parent&&Xs.parent.provides;e===a&&(a=Xs.provides=Object.create(e)),a[s]=n}}function gn(s,n,a=!1){const e=la();if(e||Ma){let t=Ma?Ma._context.provides:e?e.parent==null||e.ce?e.vnode.appContext&&e.vnode.appContext.provides:e.parent.provides:void 0;if(t&&s in t)return t[s];if(arguments.length>1)return a&&js(n)?n.call(e&&e.proxy):n}}function Tr(){return!!(la()||Ma)}const Ar={},Rr=()=>Object.create(Ar),Lr=s=>Object.getPrototypeOf(s)===Ar;function Lh(s,n,a,e=!1){const t={},l=Rr();s.propsDefaults=Object.create(null),Mr(s,n,t,l);for(const o in s.propsOptions[0])o in t||(t[o]=void 0);a?s.props=e?t:lr(t):s.type.props?s.props=t:s.props=l,s.attrs=l}function Mh(s,n,a,e){const{props:t,attrs:l,vnode:{patchFlag:o}}=s,p=ws(t),[c]=s.propsOptions;let u=!1;if((e||o>0)&&!(o&16)){if(o&8){const r=s.vnode.dynamicProps;for(let i=0;i<r.length;i++){let h=r[i];if(Mt(s.emitsOptions,h))continue;const d=n[h];if(c)if(xs(l,h))d!==l[h]&&(l[h]=d,u=!0);else{const b=Pn(h);t[b]=kl(c,p,b,d,s,!1)}else d!==l[h]&&(l[h]=d,u=!0)}}}else{Mr(s,n,t,l)&&(u=!0);let r;for(const i in p)(!n||!xs(n,i)&&((r=_a(i))===i||!xs(n,r)))&&(c?a&&(a[i]!==void 0||a[r]!==void 0)&&(t[i]=kl(c,p,i,void 0,s,!0)):delete t[i]);if(l!==p)for(const i in l)(!n||!xs(n,i))&&(delete l[i],u=!0)}u&&Yn(s.attrs,"set","")}function Mr(s,n,a,e){const[t,l]=s.propsOptions;let o=!1,p;if(n)for(let c in n){if(ce(c))continue;const u=n[c];let r;t&&xs(t,r=Pn(c))?!l||!l.includes(r)?a[r]=u:(p||(p={}))[r]=u:Mt(s.emitsOptions,c)||(!(c in e)||u!==e[c])&&(e[c]=u,o=!0)}if(l){const c=ws(a),u=p||Ss;for(let r=0;r<l.length;r++){const i=l[r];a[i]=kl(t,c,i,u[i],s,!xs(u,i))}}return o}function kl(s,n,a,e,t,l){const o=s[a];if(o!=null){const p=xs(o,"default");if(p&&e===void 0){const c=o.default;if(o.type!==Function&&!o.skipFactory&&js(c)){const{propsDefaults:u}=t;if(a in u)e=u[a];else{const r=$e(t);e=u[a]=c.call(null,n),r()}}else e=c;t.ce&&t.ce._setProp(a,e)}o[0]&&(l&&!p?e=!1:o[1]&&(e===""||e===_a(a))&&(e=!0))}return e}const Dh=new WeakMap;function Dr(s,n,a=!1){const e=a?Dh:n.propsCache,t=e.get(s);if(t)return t;const l=s.props,o={},p=[];let c=!1;if(!js(s)){const r=i=>{c=!0;const[h,d]=Dr(i,n,!0);Gs(o,h),d&&p.push(...d)};!a&&n.mixins.length&&n.mixins.forEach(r),s.extends&&r(s.extends),s.mixins&&s.mixins.forEach(r)}if(!l&&!c)return Rs(s)&&e.set(s,Ua),Ua;if(is(l))for(let r=0;r<l.length;r++){const i=Pn(l[r]);sp(i)&&(o[i]=Ss)}else if(l)for(const r in l){const i=Pn(r);if(sp(i)){const h=l[r],d=o[i]=is(h)||js(h)?{type:h}:Gs({},h),b=d.type;let f=!1,C=!0;if(is(b))for(let m=0;m<b.length;++m){const g=b[m],j=js(g)&&g.name;if(j==="Boolean"){f=!0;break}else j==="String"&&(C=!1)}else f=js(b)&&b.name==="Boolean";d[0]=f,d[1]=C,(f||xs(d,"default"))&&p.push(i)}}const u=[o,p];return Rs(s)&&e.set(s,u),u}function sp(s){return s[0]!=="$"&&!ce(s)}const fo=s=>s==="_"||s==="_ctx"||s==="$stable",bo=s=>is(s)?s.map(Fn):[Fn(s)],Oh=(s,n,a)=>{if(n._n)return n;const e=jn((...t)=>bo(n(...t)),a);return e._c=!1,e},Or=(s,n,a)=>{const e=s._ctx;for(const t in s){if(fo(t))continue;const l=s[t];if(js(l))n[t]=Oh(t,l,e);else if(l!=null){const o=bo(l);n[t]=()=>o}}},Ir=(s,n)=>{const a=bo(n);s.slots.default=()=>a},Nr=(s,n,a)=>{for(const e in n)(a||!fo(e))&&(s[e]=n[e])},Ih=(s,n,a)=>{const e=s.slots=Rr();if(s.vnode.shapeFlag&32){const t=n._;t?(Nr(e,n,a),a&&$c(e,"_",t,!0)):Or(n,e)}else n&&Ir(s,n)},Nh=(s,n,a)=>{const{vnode:e,slots:t}=s;let l=!0,o=Ss;if(e.shapeFlag&32){const p=n._;p?a&&p===1?l=!1:Nr(t,n,a):(l=!n.$stable,Or(n,t)),o=n}else n&&(Ir(s,n),o={default:1});if(l)for(const p in t)!fo(p)&&o[p]==null&&delete t[p]},pn=Zh;function Fh(s){return $h(s)}function $h(s,n){const a=At();a.__VUE__=!0;const{insert:e,remove:t,patchProp:l,createElement:o,createText:p,createComment:c,setText:u,setElementText:r,parentNode:i,nextSibling:h,setScopeId:d=$n,insertStaticContent:b}=s,f=(w,k,A,F=null,q=null,B=null,_=void 0,v=null,T=!!k.dynamicChildren)=>{if(w===k)return;w&&!Ta(w,k)&&(F=L(w),$s(w,q,B,!0),w=null),k.patchFlag===-2&&(T=!1,k.dynamicChildren=null);const{type:D,ref:Z,shapeFlag:K}=k;switch(D){case Fe:C(w,k,A,F);break;case un:m(w,k,A,F);break;case tt:w==null&&g(k,A,F,_);break;case gs:V(w,k,A,F,q,B,_,v,T);break;default:K&1?x(w,k,A,F,q,B,_,v,T):K&6?z(w,k,A,F,q,B,_,v,T):(K&64||K&128)&&D.process(w,k,A,F,q,B,_,v,T,J)}Z!=null&&q?he(Z,w&&w.ref,B,k||w,!k):Z==null&&w&&w.ref!=null&&he(w.ref,null,B,w,!0)},C=(w,k,A,F)=>{if(w==null)e(k.el=p(k.children),A,F);else{const q=k.el=w.el;k.children!==w.children&&u(q,k.children)}},m=(w,k,A,F)=>{w==null?e(k.el=c(k.children||""),A,F):k.el=w.el},g=(w,k,A,F)=>{[w.el,w.anchor]=b(w.children,k,A,F,w.el,w.anchor)},j=({el:w,anchor:k},A,F)=>{let q;for(;w&&w!==k;)q=h(w),e(w,A,F),w=q;e(k,A,F)},y=({el:w,anchor:k})=>{let A;for(;w&&w!==k;)A=h(w),t(w),w=A;t(k)},x=(w,k,A,F,q,B,_,v,T)=>{if(k.type==="svg"?_="svg":k.type==="math"&&(_="mathml"),w==null)R(k,A,F,q,B,_,v,T);else{const D=w.el&&w.el._isVueCE?w.el:null;try{D&&D._beginPatch(),P(w,k,q,B,_,v,T)}finally{D&&D._endPatch()}}},R=(w,k,A,F,q,B,_,v)=>{let T,D;const{props:Z,shapeFlag:K,transition:as,dirs:es}=w;if(T=w.el=o(w.type,B,Z&&Z.is,Z),K&8?r(T,w.children):K&16&&O(w.children,T,null,F,q,sl(w,B),_,v),es&&ka(w,null,F,"created"),E(T,w,w.scopeId,_,F),Z){for(const $ in Z)$!=="value"&&!ce($)&&l(T,$,null,Z[$],B,F);"value"in Z&&l(T,"value",null,Z.value,B),(D=Z.onVnodeBeforeMount)&&On(D,F,w)}es&&ka(w,null,F,"beforeMount");const M=Bh(q,as);M&&as.beforeEnter(T),e(T,k,A),((D=Z&&Z.onVnodeMounted)||M||es)&&pn(()=>{D&&On(D,F,w),M&&as.enter(T),es&&ka(w,null,F,"mounted")},q)},E=(w,k,A,F,q)=>{if(A&&d(w,A),F)for(let B=0;B<F.length;B++)d(w,F[B]);if(q){let B=q.subTree;if(k===B||qr(B.type)&&(B.ssContent===k||B.ssFallback===k)){const _=q.vnode;E(w,_,_.scopeId,_.slotScopeIds,q.parent)}}},O=(w,k,A,F,q,B,_,v,T=0)=>{for(let D=T;D<w.length;D++){const Z=w[D]=v?ua(w[D]):Fn(w[D]);f(null,Z,k,A,F,q,B,_,v)}},P=(w,k,A,F,q,B,_)=>{const v=k.el=w.el;let{patchFlag:T,dynamicChildren:D,dirs:Z}=k;T|=w.patchFlag&16;const K=w.props||Ss,as=k.props||Ss;let es;if(A&&xa(A,!1),(es=as.onVnodeBeforeUpdate)&&On(es,A,k,w),Z&&ka(k,w,A,"beforeUpdate"),A&&xa(A,!0),(K.innerHTML&&as.innerHTML==null||K.textContent&&as.textContent==null)&&r(v,""),D?U(w.dynamicChildren,D,v,A,F,sl(k,q),B):_||us(w,k,v,null,A,F,sl(k,q),B,!1),T>0){if(T&16)Y(v,K,as,A,q);else if(T&2&&K.class!==as.class&&l(v,"class",null,as.class,q),T&4&&l(v,"style",K.style,as.style,q),T&8){const M=k.dynamicProps;for(let $=0;$<M.length;$++){const ps=M[$],_s=K[ps],Bs=as[ps];(Bs!==_s||ps==="value")&&l(v,ps,_s,Bs,q,A)}}T&1&&w.children!==k.children&&r(v,k.children)}else!_&&D==null&&Y(v,K,as,A,q);((es=as.onVnodeUpdated)||Z)&&pn(()=>{es&&On(es,A,k,w),Z&&ka(k,w,A,"updated")},F)},U=(w,k,A,F,q,B,_)=>{for(let v=0;v<k.length;v++){const T=w[v],D=k[v],Z=T.el&&(T.type===gs||!Ta(T,D)||T.shapeFlag&198)?i(T.el):A;f(T,D,Z,null,F,q,B,_,!0)}},Y=(w,k,A,F,q)=>{if(k!==A){if(k!==Ss)for(const B in k)!ce(B)&&!(B in A)&&l(w,B,k[B],null,q,F);for(const B in A){if(ce(B))continue;const _=A[B],v=k[B];_!==v&&B!=="value"&&l(w,B,v,_,q,F)}"value"in A&&l(w,"value",k.value,A.value,q)}},V=(w,k,A,F,q,B,_,v,T)=>{const D=k.el=w?w.el:p(""),Z=k.anchor=w?w.anchor:p("");let{patchFlag:K,dynamicChildren:as,slotScopeIds:es}=k;es&&(v=v?v.concat(es):es),w==null?(e(D,A,F),e(Z,A,F),O(k.children||[],A,Z,q,B,_,v,T)):K>0&&K&64&&as&&w.dynamicChildren?(U(w.dynamicChildren,as,A,q,B,_,v),(k.key!=null||q&&k===q.subTree)&&_o(w,k,!0)):us(w,k,A,Z,q,B,_,v,T)},z=(w,k,A,F,q,B,_,v,T)=>{k.slotScopeIds=v,w==null?k.shapeFlag&512?q.ctx.activate(k,A,F,_,T):Q(k,A,F,q,B,_,T):cs(w,k,T)},Q=(w,k,A,F,q,B,_)=>{const v=w.component=od(w,F,q);if(Ie(w)&&(v.ctx.renderer=J),pd(v,!1,_),v.asyncDep){if(q&&q.registerDep(v,ls,_),!w.el){const T=v.subTree=bs(un);m(null,T,k,A),w.placeholder=T.el}}else ls(v,w,k,A,q,B,_)},cs=(w,k,A)=>{const F=k.component=w.component;if(Qh(w,k,A))if(F.asyncDep&&!F.asyncResolved){X(F,k,A);return}else F.next=k,F.update();else k.el=w.el,F.vnode=k},ls=(w,k,A,F,q,B,_)=>{const v=()=>{if(w.isMounted){let{next:K,bu:as,u:es,parent:M,vnode:$}=w;{const Ws=Fr(w);if(Ws){K&&(K.el=$.el,X(w,K,_)),Ws.asyncDep.then(()=>{w.isUnmounted||v()});return}}let ps=K,_s;xa(w,!1),K?(K.el=$.el,X(w,K,_)):K=$,as&&st(as),(_s=K.props&&K.props.onVnodeBeforeUpdate)&&On(_s,M,K,$),xa(w,!0);const Bs=ap(w),tn=w.subTree;w.subTree=Bs,f(tn,Bs,i(tn.el),L(tn),w,q,B),K.el=Bs.el,ps===null&&Jh(w,Bs.el),es&&pn(es,q),(_s=K.props&&K.props.onVnodeUpdated)&&pn(()=>On(_s,M,K,$),q)}else{let K;const{el:as,props:es}=k,{bm:M,m:$,parent:ps,root:_s,type:Bs}=w,tn=de(k);xa(w,!1),M&&st(M),!tn&&(K=es&&es.onVnodeBeforeMount)&&On(K,ps,k),xa(w,!0);{_s.ce&&_s.ce._def.shadowRoot!==!1&&_s.ce._injectChildStyle(Bs);const Ws=w.subTree=ap(w);f(null,Ws,A,F,w,q,B),k.el=Ws.el}if($&&pn($,q),!tn&&(K=es&&es.onVnodeMounted)){const Ws=k;pn(()=>On(K,ps,Ws),q)}(k.shapeFlag&256||ps&&de(ps.vnode)&&ps.vnode.shapeFlag&256)&&w.a&&pn(w.a,q),w.isMounted=!0,k=A=F=null}};w.scope.on();const T=w.effect=new Hc(v);w.scope.off();const D=w.update=T.run.bind(T),Z=w.job=T.runIfDirty.bind(T);Z.i=w,Z.id=w.uid,T.scheduler=()=>uo(Z),xa(w,!0),D()},X=(w,k,A)=>{k.component=w;const F=w.vnode.props;w.vnode=k,w.next=null,Mh(w,k.props,F,A),Nh(w,k.children,A),sa(),zo(w),na()},us=(w,k,A,F,q,B,_,v,T=!1)=>{const D=w&&w.children,Z=w?w.shapeFlag:0,K=k.children,{patchFlag:as,shapeFlag:es}=k;if(as>0){if(as&128){Qs(D,K,A,F,q,B,_,v,T);return}else if(as&256){Fs(D,K,A,F,q,B,_,v,T);return}}es&8?(Z&16&&os(D,q,B),K!==D&&r(A,K)):Z&16?es&16?Qs(D,K,A,F,q,B,_,v,T):os(D,q,B,!0):(Z&8&&r(A,""),es&16&&O(K,A,F,q,B,_,v,T))},Fs=(w,k,A,F,q,B,_,v,T)=>{w=w||Ua,k=k||Ua;const D=w.length,Z=k.length,K=Math.min(D,Z);let as;for(as=0;as<K;as++){const es=k[as]=T?ua(k[as]):Fn(k[as]);f(w[as],es,A,null,q,B,_,v,T)}D>Z?os(w,q,B,!0,!1,K):O(k,A,F,q,B,_,v,T,K)},Qs=(w,k,A,F,q,B,_,v,T)=>{let D=0;const Z=k.length;let K=w.length-1,as=Z-1;for(;D<=K&&D<=as;){const es=w[D],M=k[D]=T?ua(k[D]):Fn(k[D]);if(Ta(es,M))f(es,M,A,null,q,B,_,v,T);else break;D++}for(;D<=K&&D<=as;){const es=w[K],M=k[as]=T?ua(k[as]):Fn(k[as]);if(Ta(es,M))f(es,M,A,null,q,B,_,v,T);else break;K--,as--}if(D>K){if(D<=as){const es=as+1,M=es<Z?k[es].el:F;for(;D<=as;)f(null,k[D]=T?ua(k[D]):Fn(k[D]),A,M,q,B,_,v,T),D++}}else if(D>as)for(;D<=K;)$s(w[D],q,B,!0),D++;else{const es=D,M=D,$=new Map;for(D=M;D<=as;D++){const Js=k[D]=T?ua(k[D]):Fn(k[D]);Js.key!=null&&$.set(Js.key,D)}let ps,_s=0;const Bs=as-M+1;let tn=!1,Ws=0;const fn=new Array(Bs);for(D=0;D<Bs;D++)fn[D]=0;for(D=es;D<=K;D++){const Js=w[D];if(_s>=Bs){$s(Js,q,B,!0);continue}let Dn;if(Js.key!=null)Dn=$.get(Js.key);else for(ps=M;ps<=as;ps++)if(fn[ps-M]===0&&Ta(Js,k[ps])){Dn=ps;break}Dn===void 0?$s(Js,q,B,!0):(fn[Dn-M]=D+1,Dn>=Ws?Ws=Dn:tn=!0,f(Js,k[Dn],A,null,q,B,_,v,T),_s++)}const Ve=tn?qh(fn):Ua;for(ps=Ve.length-1,D=Bs-1;D>=0;D--){const Js=M+D,Dn=k[Js],Io=k[Js+1],No=Js+1<Z?Io.el||Io.placeholder:F;fn[D]===0?f(null,Dn,A,No,q,B,_,v,T):tn&&(ps<0||D!==Ve[ps]?Ms(Dn,A,No,2):ps--)}}},Ms=(w,k,A,F,q=null)=>{const{el:B,type:_,transition:v,children:T,shapeFlag:D}=w;if(D&6){Ms(w.component.subTree,k,A,F);return}if(D&128){w.suspense.move(k,A,F);return}if(D&64){_.move(w,k,A,J);return}if(_===gs){e(B,k,A);for(let K=0;K<T.length;K++)Ms(T[K],k,A,F);e(w.anchor,k,A);return}if(_===tt){j(w,k,A);return}if(F!==2&&D&1&&v)if(F===0)v.beforeEnter(B),e(B,k,A),pn(()=>v.enter(B),q);else{const{leave:K,delayLeave:as,afterLeave:es}=v,M=()=>{w.ctx.isUnmounted?t(B):e(B,k,A)},$=()=>{B._isLeaving&&B[Xn](!0),K(B,()=>{M(),es&&es()})};as?as(B,M,$):$()}else e(B,k,A)},$s=(w,k,A,F=!1,q=!1)=>{const{type:B,props:_,ref:v,children:T,dynamicChildren:D,shapeFlag:Z,patchFlag:K,dirs:as,cacheIndex:es}=w;if(K===-2&&(q=!1),v!=null&&(sa(),he(v,null,A,w,!0),na()),es!=null&&(k.renderCache[es]=void 0),Z&256){k.ctx.deactivate(w);return}const M=Z&1&&as,$=!de(w);let ps;if($&&(ps=_&&_.onVnodeBeforeUnmount)&&On(ps,k,w),Z&6)Mn(w.component,A,F);else{if(Z&128){w.suspense.unmount(A,F);return}M&&ka(w,null,k,"beforeUnmount"),Z&64?w.type.remove(w,k,A,J,F):D&&!D.hasOnce&&(B!==gs||K>0&&K&64)?os(D,k,A,!1,!0):(B===gs&&K&384||!q&&Z&16)&&os(T,k,A),F&&dn(w)}($&&(ps=_&&_.onVnodeUnmounted)||M)&&pn(()=>{ps&&On(ps,k,w),M&&ka(w,null,k,"unmounted")},A)},dn=w=>{const{type:k,el:A,anchor:F,transition:q}=w;if(k===gs){en(A,F);return}if(k===tt){y(w);return}const B=()=>{t(A),q&&!q.persisted&&q.afterLeave&&q.afterLeave()};if(w.shapeFlag&1&&q&&!q.persisted){const{leave:_,delayLeave:v}=q,T=()=>_(A,B);v?v(w.el,B,T):T()}else B()},en=(w,k)=>{let A;for(;w!==k;)A=h(w),t(w),w=A;t(k)},Mn=(w,k,A)=>{const{bum:F,scope:q,job:B,subTree:_,um:v,m:T,a:D}=w;np(T),np(D),F&&st(F),q.stop(),B&&(B.flags|=8,$s(_,w,k,A)),v&&pn(v,k),pn(()=>{w.isUnmounted=!0},k)},os=(w,k,A,F=!1,q=!1,B=0)=>{for(let _=B;_<w.length;_++)$s(w[_],k,A,F,q)},L=w=>{if(w.shapeFlag&6)return L(w.component.subTree);if(w.shapeFlag&128)return w.suspense.next();const k=h(w.anchor||w.el),A=k&&k[dr];return A?h(A):k};let W=!1;const H=(w,k,A)=>{w==null?k._vnode&&$s(k._vnode,null,null,!0):f(k._vnode||null,w,k,null,null,null,A),k._vnode=w,W||(W=!0,zo(),ir(),W=!1)},J={p:f,um:$s,m:Ms,r:dn,mt:Q,mc:O,pc:us,pbc:U,n:L,o:s};return{render:H,hydrate:void 0,createApp:Rh(H)}}function sl({type:s,props:n},a){return a==="svg"&&s==="foreignObject"||a==="mathml"&&s==="annotation-xml"&&n&&n.encoding&&n.encoding.includes("html")?void 0:a}function xa({effect:s,job:n},a){a?(s.flags|=32,n.flags|=4):(s.flags&=-33,n.flags&=-5)}function Bh(s,n){return(!s||s&&!s.pendingBranch)&&n&&!n.persisted}function _o(s,n,a=!1){const e=s.children,t=n.children;if(is(e)&&is(t))for(let l=0;l<e.length;l++){const o=e[l];let p=t[l];p.shapeFlag&1&&!p.dynamicChildren&&((p.patchFlag<=0||p.patchFlag===32)&&(p=t[l]=ua(t[l]),p.el=o.el),!a&&p.patchFlag!==-2&&_o(o,p)),p.type===Fe&&p.patchFlag!==-1&&(p.el=o.el),p.type===un&&!p.el&&(p.el=o.el)}}function qh(s){const n=s.slice(),a=[0];let e,t,l,o,p;const c=s.length;for(e=0;e<c;e++){const u=s[e];if(u!==0){if(t=a[a.length-1],s[t]<u){n[e]=t,a.push(e);continue}for(l=0,o=a.length-1;l<o;)p=l+o>>1,s[a[p]]<u?l=p+1:o=p;u<s[a[l]]&&(l>0&&(n[e]=a[l-1]),a[l]=e)}}for(l=a.length,o=a[l-1];l-- >0;)a[l]=o,o=n[o];return a}function Fr(s){const n=s.subTree.component;if(n)return n.asyncDep&&!n.asyncResolved?n:Fr(n)}function np(s){if(s)for(let n=0;n<s.length;n++)s[n].flags|=8}const zh=Symbol.for("v-scx"),Vh=()=>gn(zh);function Uh(s,n){return yo(s,null,n)}function Vs(s,n,a){return yo(s,n,a)}function yo(s,n,a=Ss){const{immediate:e,deep:t,flush:l,once:o}=a,p=Gs({},a),c=n&&e||!n&&l!=="post";let u;if(Ka){if(l==="sync"){const d=Vh();u=d.__watcherHandles||(d.__watcherHandles=[])}else if(!c){const d=()=>{};return d.stop=$n,d.resume=$n,d.pause=$n,d}}const r=Xs;p.call=(d,b,f)=>Rn(d,r,b,f);let i=!1;l==="post"?p.scheduler=d=>{pn(d,r&&r.suspense)}:l!=="sync"&&(i=!0,p.scheduler=(d,b)=>{b?d():uo(d)}),p.augmentJob=d=>{n&&(d.flags|=4),i&&(d.flags|=2,r&&(d.id=r.uid,d.i=r))};const h=th(s,n,p);return Ka&&(u?u.push(h):c&&h()),h}function Hh(s,n,a){const e=this.proxy,t=Is(s)?s.includes(".")?$r(e,s):()=>e[s]:s.bind(e,e);let l;js(n)?l=n:(l=n.handler,a=n);const o=$e(this),p=yo(t,l.bind(e),a);return o(),p}function $r(s,n){const a=n.split(".");return()=>{let e=s;for(let t=0;t<a.length&&e;t++)e=e[a[t]];return e}}const Gh=(s,n)=>n==="modelValue"||n==="model-value"?s.modelModifiers:s[`${n}Modifiers`]||s[`${Pn(n)}Modifiers`]||s[`${_a(n)}Modifiers`];function Wh(s,n,...a){if(s.isUnmounted)return;const e=s.vnode.props||Ss;let t=a;const l=n.startsWith("update:"),o=l&&Gh(e,n.slice(7));o&&(o.trim&&(t=a.map(r=>Is(r)?r.trim():r)),o.number&&(t=a.map(so)));let p,c=e[p=Wt(n)]||e[p=Wt(Pn(n))];!c&&l&&(c=e[p=Wt(_a(n))]),c&&Rn(c,s,6,t);const u=e[p+"Once"];if(u){if(!s.emitted)s.emitted={};else if(s.emitted[p])return;s.emitted[p]=!0,Rn(u,s,6,t)}}const Kh=new WeakMap;function Br(s,n,a=!1){const e=a?Kh:n.emitsCache,t=e.get(s);if(t!==void 0)return t;const l=s.emits;let o={},p=!1;if(!js(s)){const c=u=>{const r=Br(u,n,!0);r&&(p=!0,Gs(o,r))};!a&&n.mixins.length&&n.mixins.forEach(c),s.extends&&c(s.extends),s.mixins&&s.mixins.forEach(c)}return!l&&!p?(Rs(s)&&e.set(s,null),null):(is(l)?l.forEach(c=>o[c]=null):Gs(o,l),Rs(s)&&e.set(s,o),o)}function Mt(s,n){return!s||!St(n)?!1:(n=n.slice(2).replace(/Once$/,""),xs(s,n[0].toLowerCase()+n.slice(1))||xs(s,_a(n))||xs(s,n))}function ap(s){const{type:n,vnode:a,proxy:e,withProxy:t,propsOptions:[l],slots:o,attrs:p,emit:c,render:u,renderCache:r,props:i,data:h,setupState:d,ctx:b,inheritAttrs:f}=s,C=dt(s);let m,g;try{if(a.shapeFlag&4){const y=t||e,x=y;m=Fn(u.call(x,y,r,i,d,h,b)),g=p}else{const y=n;m=Fn(y.length>1?y(i,{attrs:p,slots:o,emit:c}):y(i,null)),g=n.props?p:Xh(p)}}catch(y){je.length=0,Oe(y,s,1),m=bs(un)}let j=m;if(g&&f!==!1){const y=Object.keys(g),{shapeFlag:x}=j;y.length&&x&7&&(l&&y.some(Ql)&&(g=Yh(g,l)),j=ga(j,g,!1,!0))}return a.dirs&&(j=ga(j,null,!1,!0),j.dirs=j.dirs?j.dirs.concat(a.dirs):a.dirs),a.transition&&Ce(j,a.transition),m=j,dt(C),m}const Xh=s=>{let n;for(const a in s)(a==="class"||a==="style"||St(a))&&((n||(n={}))[a]=s[a]);return n},Yh=(s,n)=>{const a={};for(const e in s)(!Ql(e)||!(e.slice(9)in n))&&(a[e]=s[e]);return a};function Qh(s,n,a){const{props:e,children:t,component:l}=s,{props:o,children:p,patchFlag:c}=n,u=l.emitsOptions;if(n.dirs||n.transition)return!0;if(a&&c>=0){if(c&1024)return!0;if(c&16)return e?ep(e,o,u):!!o;if(c&8){const r=n.dynamicProps;for(let i=0;i<r.length;i++){const h=r[i];if(o[h]!==e[h]&&!Mt(u,h))return!0}}}else return(t||p)&&(!p||!p.$stable)?!0:e===o?!1:e?o?ep(e,o,u):!0:!!o;return!1}function ep(s,n,a){const e=Object.keys(n);if(e.length!==Object.keys(s).length)return!0;for(let t=0;t<e.length;t++){const l=e[t];if(n[l]!==s[l]&&!Mt(a,l))return!0}return!1}function Jh({vnode:s,parent:n},a){for(;n;){const e=n.subTree;if(e.suspense&&e.suspense.activeBranch===s&&(e.el=s.el),e===s)(s=n.vnode).el=a,n=n.parent;else break}}const qr=s=>s.__isSuspense;function Zh(s,n){n&&n.pendingBranch?is(s)?n.effects.push(...s):n.effects.push(s):ph(s)}const gs=Symbol.for("v-fgt"),Fe=Symbol.for("v-txt"),un=Symbol.for("v-cmt"),tt=Symbol.for("v-stc"),je=[];let yn=null;function I(s=!1){je.push(yn=s?null:[])}function sd(){je.pop(),yn=je[je.length-1]||null}let ke=1;function gt(s,n=!1){ke+=s,s<0&&yn&&n&&(yn.hasOnce=!0)}function zr(s){return s.dynamicChildren=ke>0?yn||Ua:null,sd(),ke>0&&yn&&yn.push(s),s}function N(s,n,a,e,t,l){return zr(S(s,n,a,e,t,l,!0))}function Aa(s,n,a,e,t){return zr(bs(s,n,a,e,t,!0))}function ft(s){return s?s.__v_isVNode===!0:!1}function Ta(s,n){return s.type===n.type&&s.key===n.key}const Vr=({key:s})=>s??null,lt=({ref:s,ref_key:n,ref_for:a})=>(typeof s=="number"&&(s=""+s),s!=null?Is(s)||Os(s)||js(s)?{i:_n,r:s,k:n,f:!!a}:s:null);function S(s,n=null,a=null,e=0,t=null,l=s===gs?0:1,o=!1,p=!1){const c={__v_isVNode:!0,__v_skip:!0,type:s,props:n,key:n&&Vr(n),ref:n&&lt(n),scopeId:hr,slotScopeIds:null,children:a,component:null,suspense:null,ssContent:null,ssFallback:null,dirs:null,transition:null,el:null,anchor:null,target:null,targetStart:null,targetAnchor:null,staticCount:0,shapeFlag:l,patchFlag:e,dynamicProps:t,dynamicChildren:null,appContext:null,ctx:_n};return p?(vo(c,a),l&128&&s.normalize(c)):a&&(c.shapeFlag|=Is(a)?8:16),ke>0&&!o&&yn&&(c.patchFlag>0||l&6)&&c.patchFlag!==32&&yn.push(c),c}const bs=nd;function nd(s,n=null,a=null,e=0,t=null,l=!1){if((!s||s===xr)&&(s=un),ft(s)){const p=ga(s,n,!0);return a&&vo(p,a),ke>0&&!l&&yn&&(p.shapeFlag&6?yn[yn.indexOf(s)]=p:yn.push(p)),p.patchFlag=-2,p}if(hd(s)&&(s=s.__vccOpts),n){n=ad(n);let{class:p,style:c}=n;p&&!Is(p)&&(n.class=Ts(p)),Rs(c)&&(co(c)&&!is(c)&&(c=Gs({},c)),n.style=qn(c))}const o=Is(s)?1:qr(s)?128:mr(s)?64:Rs(s)?4:js(s)?2:0;return S(s,n,a,e,t,o,l,!0)}function ad(s){return s?co(s)||Lr(s)?Gs({},s):s:null}function ga(s,n,a=!1,e=!1){const{props:t,ref:l,patchFlag:o,children:p,transition:c}=s,u=n?ed(t||{},n):t,r={__v_isVNode:!0,__v_skip:!0,type:s.type,props:u,key:u&&Vr(u),ref:n&&n.ref?a&&l?is(l)?l.concat(lt(n)):[l,lt(n)]:lt(n):l,scopeId:s.scopeId,slotScopeIds:s.slotScopeIds,children:p,target:s.target,targetStart:s.targetStart,targetAnchor:s.targetAnchor,staticCount:s.staticCount,shapeFlag:s.shapeFlag,patchFlag:n&&s.type!==gs?o===-1?16:o|16:o,dynamicProps:s.dynamicProps,dynamicChildren:s.dynamicChildren,appContext:s.appContext,dirs:s.dirs,transition:c,component:s.component,suspense:s.suspense,ssContent:s.ssContent&&ga(s.ssContent),ssFallback:s.ssFallback&&ga(s.ssFallback),placeholder:s.placeholder,el:s.el,anchor:s.anchor,ctx:s.ctx,ce:s.ce};return c&&e&&Ce(r,c.clone(r)),r}function Cs(s=" ",n=0){return bs(Fe,null,s,n)}function tp(s,n){const a=bs(tt,null,s);return a.staticCount=n,a}function ms(s="",n=!1){return n?(I(),Aa(un,null,s)):bs(un,null,s)}function Fn(s){return s==null||typeof s=="boolean"?bs(un):is(s)?bs(gs,null,s.slice()):ft(s)?ua(s):bs(Fe,null,String(s))}function ua(s){return s.el===null&&s.patchFlag!==-1||s.memo?s:ga(s)}function vo(s,n){let a=0;const{shapeFlag:e}=s;if(n==null)n=null;else if(is(n))a=16;else if(typeof n=="object")if(e&65){const t=n.default;t&&(t._c&&(t._d=!1),vo(s,t()),t._c&&(t._d=!0));return}else{a=32;const t=n._;!t&&!Lr(n)?n._ctx=_n:t===3&&_n&&(_n.slots._===1?n._=1:(n._=2,s.patchFlag|=1024))}else js(n)?(n={default:n,_ctx:_n},a=32):(n=String(n),e&64?(a=16,n=[Cs(n)]):a=8);s.children=n,s.shapeFlag|=a}function ed(...s){const n={};for(let a=0;a<s.length;a++){const e=s[a];for(const t in e)if(t==="class")n.class!==e.class&&(n.class=Ts([n.class,e.class]));else if(t==="style")n.style=qn([n.style,e.style]);else if(St(t)){const l=n[t],o=e[t];o&&l!==o&&!(is(l)&&l.includes(o))&&(n[t]=l?[].concat(l,o):o)}else t!==""&&(n[t]=e[t])}return n}function On(s,n,a,e=null){Rn(s,n,7,[a,e])}const td=Er();let ld=0;function od(s,n,a){const e=s.type,t=(n?n.appContext:s.appContext)||td,l={uid:ld++,vnode:s,type:e,parent:n,appContext:t,root:null,next:null,subTree:null,effect:null,update:null,job:null,scope:new Vc(!0),render:null,proxy:null,exposed:null,exposeProxy:null,withProxy:null,provides:n?n.provides:Object.create(t.provides),ids:n?n.ids:["",0,0],accessCache:null,renderCache:[],components:null,directives:null,propsOptions:Dr(e,t),emitsOptions:Br(e,t),emit:null,emitted:null,propsDefaults:Ss,inheritAttrs:e.inheritAttrs,ctx:Ss,data:Ss,props:Ss,attrs:Ss,slots:Ss,refs:Ss,setupState:Ss,setupContext:null,suspense:a,suspenseId:a?a.pendingId:0,asyncDep:null,asyncResolved:!1,isMounted:!1,isUnmounted:!1,isDeactivated:!1,bc:null,c:null,bm:null,m:null,bu:null,u:null,um:null,bum:null,da:null,a:null,rtg:null,rtc:null,ec:null,sp:null};return l.ctx={_:l},l.root=n?n.root:l,l.emit=Wh.bind(null,l),s.ce&&s.ce(l),l}let Xs=null;const la=()=>Xs||_n;let bt,xl;{const s=At(),n=(a,e)=>{let t;return(t=s[a])||(t=s[a]=[]),t.push(e),l=>{t.length>1?t.forEach(o=>o(l)):t[0](l)}};bt=n("__VUE_INSTANCE_SETTERS__",a=>Xs=a),xl=n("__VUE_SSR_SETTERS__",a=>Ka=a)}const $e=s=>{const n=Xs;return bt(s),s.scope.on(),()=>{s.scope.off(),bt(n)}},lp=()=>{Xs&&Xs.scope.off(),bt(null)};function Ur(s){return s.vnode.shapeFlag&4}let Ka=!1;function pd(s,n=!1,a=!1){n&&xl(n);const{props:e,children:t}=s.vnode,l=Ur(s);Lh(s,e,l,n),Ih(s,t,a||n);const o=l?cd(s,n):void 0;return n&&xl(!1),o}function cd(s,n){const a=s.type;s.accessCache=Object.create(null),s.proxy=new Proxy(s.ctx,kh);const{setup:e}=a;if(e){sa();const t=s.setupContext=e.length>1?id(s):null,l=$e(s),o=De(e,s,0,[s.props,t]),p=Ic(o);if(na(),l(),(p||s.sp)&&!de(s)&&ho(s),p){if(o.then(lp,lp),n)return o.then(c=>{op(s,c)}).catch(c=>{Oe(c,s,0)});s.asyncDep=o}else op(s,o)}else Hr(s)}function op(s,n,a){js(n)?s.type.__ssrInlineRender?s.ssrRender=n:s.render=n:Rs(n)&&(s.setupState=pr(n)),Hr(s)}function Hr(s,n,a){const e=s.type;s.render||(s.render=e.render||$n);{const t=$e(s);sa();try{xh(s)}finally{na(),t()}}}const rd={get(s,n){return nn(s,"get",""),s[n]}};function id(s){const n=a=>{s.exposed=a||{}};return{attrs:new Proxy(s.attrs,rd),slots:s.slots,emit:s.emit,expose:n}}function Dt(s){return s.exposed?s.exposeProxy||(s.exposeProxy=new Proxy(pr(ro(s.exposed)),{get(n,a){if(a in n)return n[a];if(a in me)return me[a](s)},has(n,a){return a in n||a in me}})):s.proxy}function ud(s,n=!0){return js(s)?s.displayName||s.name:s.name||n&&s.__name}function hd(s){return js(s)&&"__vccOpts"in s}const ns=(s,n)=>ah(s,n,Ka);function Be(s,n,a){try{gt(-1);const e=arguments.length;return e===2?Rs(n)&&!is(n)?ft(n)?bs(s,null,[n]):bs(s,n):bs(s,null,n):(e>3?a=Array.prototype.slice.call(arguments,2):e===3&&ft(a)&&(a=[a]),bs(s,n,a))}finally{gt(1)}}const dd="3.5.24";/**
* @vue/runtime-dom v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/let Sl;const pp=typeof window<"u"&&window.trustedTypes;if(pp)try{Sl=pp.createPolicy("vue",{createHTML:s=>s})}catch{}const Gr=Sl?s=>Sl.createHTML(s):s=>s,md="http://www.w3.org/2000/svg",jd="http://www.w3.org/1998/Math/MathML",Kn=typeof document<"u"?document:null,cp=Kn&&Kn.createElement("template"),gd={insert:(s,n,a)=>{n.insertBefore(s,a||null)},remove:s=>{const n=s.parentNode;n&&n.removeChild(s)},createElement:(s,n,a,e)=>{const t=n==="svg"?Kn.createElementNS(md,s):n==="mathml"?Kn.createElementNS(jd,s):a?Kn.createElement(s,{is:a}):Kn.createElement(s);return s==="select"&&e&&e.multiple!=null&&t.setAttribute("multiple",e.multiple),t},createText:s=>Kn.createTextNode(s),createComment:s=>Kn.createComment(s),setText:(s,n)=>{s.nodeValue=n},setElementText:(s,n)=>{s.textContent=n},parentNode:s=>s.parentNode,nextSibling:s=>s.nextSibling,querySelector:s=>Kn.querySelector(s),setScopeId(s,n){s.setAttribute(n,"")},insertStaticContent(s,n,a,e,t,l){const o=a?a.previousSibling:n.lastChild;if(t&&(t===l||t.nextSibling))for(;n.insertBefore(t.cloneNode(!0),a),!(t===l||!(t=t.nextSibling)););else{cp.innerHTML=Gr(e==="svg"?`<svg>${s}</svg>`:e==="mathml"?`<math>${s}</math>`:s);const p=cp.content;if(e==="svg"||e==="mathml"){const c=p.firstChild;for(;c.firstChild;)p.appendChild(c.firstChild);p.removeChild(c)}n.insertBefore(p,a)}return[o?o.nextSibling:n.firstChild,a?a.previousSibling:n.lastChild]}},oa="transition",ae="animation",xe=Symbol("_vtc"),Wr={name:String,type:String,css:{type:Boolean,default:!0},duration:[String,Number,Object],enterFromClass:String,enterActiveClass:String,enterToClass:String,appearFromClass:String,appearActiveClass:String,appearToClass:String,leaveFromClass:String,leaveActiveClass:String,leaveToClass:String},fd=Gs({},fr,Wr),bd=s=>(s.displayName="Transition",s.props=fd,s),_d=bd((s,{slots:n})=>Be(uh,yd(s),n)),Sa=(s,n=[])=>{is(s)?s.forEach(a=>a(...n)):s&&s(...n)},rp=s=>s?is(s)?s.some(n=>n.length>1):s.length>1:!1;function yd(s){const n={};for(const V in s)V in Wr||(n[V]=s[V]);if(s.css===!1)return n;const{name:a="v",type:e,duration:t,enterFromClass:l=`${a}-enter-from`,enterActiveClass:o=`${a}-enter-active`,enterToClass:p=`${a}-enter-to`,appearFromClass:c=l,appearActiveClass:u=o,appearToClass:r=p,leaveFromClass:i=`${a}-leave-from`,leaveActiveClass:h=`${a}-leave-active`,leaveToClass:d=`${a}-leave-to`}=s,b=vd(t),f=b&&b[0],C=b&&b[1],{onBeforeEnter:m,onEnter:g,onEnterCancelled:j,onLeave:y,onLeaveCancelled:x,onBeforeAppear:R=m,onAppear:E=g,onAppearCancelled:O=j}=n,P=(V,z,Q,cs)=>{V._enterCancelled=cs,Pa(V,z?r:p),Pa(V,z?u:o),Q&&Q()},U=(V,z)=>{V._isLeaving=!1,Pa(V,i),Pa(V,d),Pa(V,h),z&&z()},Y=V=>(z,Q)=>{const cs=V?E:g,ls=()=>P(z,V,Q);Sa(cs,[z,ls]),ip(()=>{Pa(z,V?c:l),Hn(z,V?r:p),rp(cs)||up(z,e,f,ls)})};return Gs(n,{onBeforeEnter(V){Sa(m,[V]),Hn(V,l),Hn(V,o)},onBeforeAppear(V){Sa(R,[V]),Hn(V,c),Hn(V,u)},onEnter:Y(!1),onAppear:Y(!0),onLeave(V,z){V._isLeaving=!0;const Q=()=>U(V,z);Hn(V,i),V._enterCancelled?(Hn(V,h),mp(V)):(mp(V),Hn(V,h)),ip(()=>{V._isLeaving&&(Pa(V,i),Hn(V,d),rp(y)||up(V,e,C,Q))}),Sa(y,[V,Q])},onEnterCancelled(V){P(V,!1,void 0,!0),Sa(j,[V])},onAppearCancelled(V){P(V,!0,void 0,!0),Sa(O,[V])},onLeaveCancelled(V){U(V),Sa(x,[V])}})}function vd(s){if(s==null)return null;if(Rs(s))return[nl(s.enter),nl(s.leave)];{const n=nl(s);return[n,n]}}function nl(s){return wu(s)}function Hn(s,n){n.split(/\s+/).forEach(a=>a&&s.classList.add(a)),(s[xe]||(s[xe]=new Set)).add(n)}function Pa(s,n){n.split(/\s+/).forEach(e=>e&&s.classList.remove(e));const a=s[xe];a&&(a.delete(n),a.size||(s[xe]=void 0))}function ip(s){requestAnimationFrame(()=>{requestAnimationFrame(s)})}let wd=0;function up(s,n,a,e){const t=s._endId=++wd,l=()=>{t===s._endId&&e()};if(a!=null)return setTimeout(l,a);const{type:o,timeout:p,propCount:c}=Cd(s,n);if(!o)return e();const u=o+"end";let r=0;const i=()=>{s.removeEventListener(u,h),l()},h=d=>{d.target===s&&++r>=c&&i()};setTimeout(()=>{r<c&&i()},p+1),s.addEventListener(u,h)}function Cd(s,n){const a=window.getComputedStyle(s),e=b=>(a[b]||"").split(", "),t=e(`${oa}Delay`),l=e(`${oa}Duration`),o=hp(t,l),p=e(`${ae}Delay`),c=e(`${ae}Duration`),u=hp(p,c);let r=null,i=0,h=0;n===oa?o>0&&(r=oa,i=o,h=l.length):n===ae?u>0&&(r=ae,i=u,h=c.length):(i=Math.max(o,u),r=i>0?o>u?oa:ae:null,h=r?r===oa?l.length:c.length:0);const d=r===oa&&/\b(?:transform|all)(?:,|$)/.test(e(`${oa}Property`).toString());return{type:r,timeout:i,propCount:h,hasTransform:d}}function hp(s,n){for(;s.length<n.length;)s=s.concat(s);return Math.max(...n.map((a,e)=>dp(a)+dp(s[e])))}function dp(s){return s==="auto"?0:Number(s.slice(0,-1).replace(",","."))*1e3}function mp(s){return(s?s.ownerDocument:document).body.offsetHeight}function kd(s,n,a){const e=s[xe];e&&(n=(n?[n,...e]:[...e]).join(" ")),n==null?s.removeAttribute("class"):a?s.setAttribute("class",n):s.className=n}const _t=Symbol("_vod"),Kr=Symbol("_vsh"),yt={name:"show",beforeMount(s,{value:n},{transition:a}){s[_t]=s.style.display==="none"?"":s.style.display,a&&n?a.beforeEnter(s):ee(s,n)},mounted(s,{value:n},{transition:a}){a&&n&&a.enter(s)},updated(s,{value:n,oldValue:a},{transition:e}){!n!=!a&&(e?n?(e.beforeEnter(s),ee(s,!0),e.enter(s)):e.leave(s,()=>{ee(s,!1)}):ee(s,n))},beforeUnmount(s,{value:n}){ee(s,n)}};function ee(s,n){s.style.display=n?s[_t]:"none",s[Kr]=!n}const xd=Symbol(""),Sd=/(?:^|;)\s*display\s*:/;function Pd(s,n,a){const e=s.style,t=Is(a);let l=!1;if(a&&!t){if(n)if(Is(n))for(const o of n.split(";")){const p=o.slice(0,o.indexOf(":")).trim();a[p]==null&&ot(e,p,"")}else for(const o in n)a[o]==null&&ot(e,o,"");for(const o in a)o==="display"&&(l=!0),ot(e,o,a[o])}else if(t){if(n!==a){const o=e[xd];o&&(a+=";"+o),e.cssText=a,l=Sd.test(a)}}else n&&s.removeAttribute("style");_t in s&&(s[_t]=l?e.display:"",s[Kr]&&(e.display="none"))}const jp=/\s*!important$/;function ot(s,n,a){if(is(a))a.forEach(e=>ot(s,n,e));else if(a==null&&(a=""),n.startsWith("--"))s.setProperty(n,a);else{const e=Ed(s,n);jp.test(a)?s.setProperty(_a(e),a.replace(jp,""),"important"):s[e]=a}}const gp=["Webkit","Moz","ms"],al={};function Ed(s,n){const a=al[n];if(a)return a;let e=Pn(n);if(e!=="filter"&&e in s)return al[n]=e;e=Tt(e);for(let t=0;t<gp.length;t++){const l=gp[t]+e;if(l in s)return al[n]=l}return n}const fp="http://www.w3.org/1999/xlink";function bp(s,n,a,e,t,l=Eu(n)){e&&n.startsWith("xlink:")?a==null?s.removeAttributeNS(fp,n.slice(6,n.length)):s.setAttributeNS(fp,n,a):a==null||l&&!Bc(a)?s.removeAttribute(n):s.setAttribute(n,l?"":ba(a)?String(a):a)}function _p(s,n,a,e,t){if(n==="innerHTML"||n==="textContent"){a!=null&&(s[n]=n==="innerHTML"?Gr(a):a);return}const l=s.tagName;if(n==="value"&&l!=="PROGRESS"&&!l.includes("-")){const p=l==="OPTION"?s.getAttribute("value")||"":s.value,c=a==null?s.type==="checkbox"?"on":"":String(a);(p!==c||!("_value"in s))&&(s.value=c),a==null&&s.removeAttribute(n),s._value=a;return}let o=!1;if(a===""||a==null){const p=typeof s[n];p==="boolean"?a=Bc(a):a==null&&p==="string"?(a="",o=!0):p==="number"&&(a=0,o=!0)}try{s[n]=a}catch{}o&&s.removeAttribute(t||n)}function qa(s,n,a,e){s.addEventListener(n,a,e)}function Td(s,n,a,e){s.removeEventListener(n,a,e)}const yp=Symbol("_vei");function Ad(s,n,a,e,t=null){const l=s[yp]||(s[yp]={}),o=l[n];if(e&&o)o.value=e;else{const[p,c]=Rd(n);if(e){const u=l[n]=Dd(e,t);qa(s,p,u,c)}else o&&(Td(s,p,o,c),l[n]=void 0)}}const vp=/(?:Once|Passive|Capture)$/;function Rd(s){let n;if(vp.test(s)){n={};let e;for(;e=s.match(vp);)s=s.slice(0,s.length-e[0].length),n[e[0].toLowerCase()]=!0}return[s[2]===":"?s.slice(3):_a(s.slice(2)),n]}let el=0;const Ld=Promise.resolve(),Md=()=>el||(Ld.then(()=>el=0),el=Date.now());function Dd(s,n){const a=e=>{if(!e._vts)e._vts=Date.now();else if(e._vts<=a.attached)return;Rn(Od(e,a.value),n,5,[e])};return a.value=s,a.attached=Md(),a}function Od(s,n){if(is(n)){const a=s.stopImmediatePropagation;return s.stopImmediatePropagation=()=>{a.call(s),s._stopped=!0},n.map(e=>t=>!t._stopped&&e&&e(t))}else return n}const wp=s=>s.charCodeAt(0)===111&&s.charCodeAt(1)===110&&s.charCodeAt(2)>96&&s.charCodeAt(2)<123,Id=(s,n,a,e,t,l)=>{const o=t==="svg";n==="class"?kd(s,e,o):n==="style"?Pd(s,a,e):St(n)?Ql(n)||Ad(s,n,a,e,l):(n[0]==="."?(n=n.slice(1),!0):n[0]==="^"?(n=n.slice(1),!1):Nd(s,n,e,o))?(_p(s,n,e),!s.tagName.includes("-")&&(n==="value"||n==="checked"||n==="selected")&&bp(s,n,e,o,l,n!=="value")):s._isVueCE&&(/[A-Z]/.test(n)||!Is(e))?_p(s,Pn(n),e,l,n):(n==="true-value"?s._trueValue=e:n==="false-value"&&(s._falseValue=e),bp(s,n,e,o))};function Nd(s,n,a,e){if(e)return!!(n==="innerHTML"||n==="textContent"||n in s&&wp(n)&&js(a));if(n==="spellcheck"||n==="draggable"||n==="translate"||n==="autocorrect"||n==="sandbox"&&s.tagName==="IFRAME"||n==="form"||n==="list"&&s.tagName==="INPUT"||n==="type"&&s.tagName==="TEXTAREA")return!1;if(n==="width"||n==="height"){const t=s.tagName;if(t==="IMG"||t==="VIDEO"||t==="CANVAS"||t==="SOURCE")return!1}return wp(n)&&Is(a)?!1:n in s}const Cp=s=>{const n=s.props["onUpdate:modelValue"]||!1;return is(n)?a=>st(n,a):n};function Fd(s){s.target.composing=!0}function kp(s){const n=s.target;n.composing&&(n.composing=!1,n.dispatchEvent(new Event("input")))}const tl=Symbol("_assign");function xp(s,n,a){return n&&(s=s.trim()),a&&(s=so(s)),s}const U1={created(s,{modifiers:{lazy:n,trim:a,number:e}},t){s[tl]=Cp(t);const l=e||t.props&&t.props.type==="number";qa(s,n?"change":"input",o=>{o.target.composing||s[tl](xp(s.value,a,l))}),(a||l)&&qa(s,"change",()=>{s.value=xp(s.value,a,l)}),n||(qa(s,"compositionstart",Fd),qa(s,"compositionend",kp),qa(s,"change",kp))},mounted(s,{value:n}){s.value=n??""},beforeUpdate(s,{value:n,oldValue:a,modifiers:{lazy:e,trim:t,number:l}},o){if(s[tl]=Cp(o),s.composing)return;const p=(l||s.type==="number")&&!/^0\d/.test(s.value)?so(s.value):s.value,c=n??"";p!==c&&(document.activeElement===s&&s.type!=="range"&&(e&&n===a||t&&s.value.trim()===c)||(s.value=c))}},$d=["ctrl","shift","alt","meta"],Bd={stop:s=>s.stopPropagation(),prevent:s=>s.preventDefault(),self:s=>s.target!==s.currentTarget,ctrl:s=>!s.ctrlKey,shift:s=>!s.shiftKey,alt:s=>!s.altKey,meta:s=>!s.metaKey,left:s=>"button"in s&&s.button!==0,middle:s=>"button"in s&&s.button!==1,right:s=>"button"in s&&s.button!==2,exact:(s,n)=>$d.some(a=>s[`${a}Key`]&&!n.includes(a))},vn=(s,n)=>{const a=s._withMods||(s._withMods={}),e=n.join(".");return a[e]||(a[e]=((t,...l)=>{for(let o=0;o<n.length;o++){const p=Bd[n[o]];if(p&&p(t,n))return}return s(t,...l)}))},qd={esc:"escape",space:" ",up:"arrow-up",left:"arrow-left",right:"arrow-right",down:"arrow-down",delete:"backspace"},Sp=(s,n)=>{const a=s._withKeys||(s._withKeys={}),e=n.join(".");return a[e]||(a[e]=(t=>{if(!("key"in t))return;const l=_a(t.key);if(n.some(o=>o===l||qd[o]===l))return s(t)}))},zd=Gs({patchProp:Id},gd);let Pp;function Vd(){return Pp||(Pp=Fh(zd))}const Ud=((...s)=>{const n=Vd().createApp(...s),{mount:a}=n;return n.mount=e=>{const t=Gd(e);if(!t)return;const l=n._component;!js(l)&&!l.render&&!l.template&&(l.template=t.innerHTML),t.nodeType===1&&(t.textContent="");const o=a(t,!1,Hd(t));return t instanceof Element&&(t.removeAttribute("v-cloak"),t.setAttribute("data-v-app","")),o},n});function Hd(s){if(s instanceof SVGElement)return"svg";if(typeof MathMLElement=="function"&&s instanceof MathMLElement)return"mathml"}function Gd(s){return Is(s)?document.querySelector(s):s}/*!
 * pinia v3.0.3
 * (c) 2025 Eduardo San Martin Morote
 * @license MIT
 */let Xr;const Ot=s=>Xr=s,Yr=Symbol();function Pl(s){return s&&typeof s=="object"&&Object.prototype.toString.call(s)==="[object Object]"&&typeof s.toJSON!="function"}var ge;(function(s){s.direct="direct",s.patchObject="patch object",s.patchFunction="patch function"})(ge||(ge={}));function Wd(){const s=no(!0),n=s.run(()=>rs({}));let a=[],e=[];const t=ro({install(l){Ot(t),t._a=l,l.provide(Yr,t),l.config.globalProperties.$pinia=t,e.forEach(o=>a.push(o)),e=[]},use(l){return this._a?a.push(l):e.push(l),this},_p:a,_a:null,_e:s,_s:new Map,state:n});return t}const Qr=()=>{};function Ep(s,n,a,e=Qr){s.push(n);const t=()=>{const l=s.indexOf(n);l>-1&&(s.splice(l,1),e())};return!a&&Uc()&&Tu(t),t}function Na(s,...n){s.slice().forEach(a=>{a(...n)})}const Kd=s=>s(),Tp=Symbol(),ll=Symbol();function El(s,n){s instanceof Map&&n instanceof Map?n.forEach((a,e)=>s.set(e,a)):s instanceof Set&&n instanceof Set&&n.forEach(s.add,s);for(const a in n){if(!n.hasOwnProperty(a))continue;const e=n[a],t=s[a];Pl(t)&&Pl(e)&&s.hasOwnProperty(a)&&!Os(e)&&!ma(e)?s[a]=El(t,e):s[a]=e}return s}const Xd=Symbol();function Yd(s){return!Pl(s)||!Object.prototype.hasOwnProperty.call(s,Xd)}const{assign:ra}=Object;function Qd(s){return!!(Os(s)&&s.effect)}function Jd(s,n,a,e){const{state:t,actions:l,getters:o}=n,p=a.state.value[s];let c;function u(){p||(a.state.value[s]=t?t():{});const r=Ju(a.state.value[s]);return ra(r,l,Object.keys(o||{}).reduce((i,h)=>(i[h]=ro(ns(()=>{Ot(a);const d=a._s.get(s);return o[h].call(d,d)})),i),{}))}return c=Jr(s,u,n,a,e,!0),c}function Jr(s,n,a={},e,t,l){let o;const p=ra({actions:{}},a),c={deep:!0};let u,r,i=[],h=[],d;const b=e.state.value[s];!l&&!b&&(e.state.value[s]={}),rs({});let f;function C(O){let P;u=r=!1,typeof O=="function"?(O(e.state.value[s]),P={type:ge.patchFunction,storeId:s,events:d}):(El(e.state.value[s],O),P={type:ge.patchObject,payload:O,storeId:s,events:d});const U=f=Symbol();bn().then(()=>{f===U&&(u=!0)}),r=!0,Na(i,P,e.state.value[s])}const m=l?function(){const{state:P}=a,U=P?P():{};this.$patch(Y=>{ra(Y,U)})}:Qr;function g(){o.stop(),i=[],h=[],e._s.delete(s)}const j=(O,P="")=>{if(Tp in O)return O[ll]=P,O;const U=function(){Ot(e);const Y=Array.from(arguments),V=[],z=[];function Q(X){V.push(X)}function cs(X){z.push(X)}Na(h,{args:Y,name:U[ll],store:x,after:Q,onError:cs});let ls;try{ls=O.apply(this&&this.$id===s?this:x,Y)}catch(X){throw Na(z,X),X}return ls instanceof Promise?ls.then(X=>(Na(V,X),X)).catch(X=>(Na(z,X),Promise.reject(X))):(Na(V,ls),ls)};return U[Tp]=!0,U[ll]=P,U},y={_p:e,$id:s,$onAction:Ep.bind(null,h),$patch:C,$reset:m,$subscribe(O,P={}){const U=Ep(i,O,P.detached,()=>Y()),Y=o.run(()=>Vs(()=>e.state.value[s],V=>{(P.flush==="sync"?r:u)&&O({storeId:s,type:ge.direct,events:d},V)},ra({},c,P)));return U},$dispose:g},x=Me(y);e._s.set(s,x);const E=(e._a&&e._a.runWithContext||Kd)(()=>e._e.run(()=>(o=no()).run(()=>n({action:j}))));for(const O in E){const P=E[O];if(Os(P)&&!Qd(P)||ma(P))l||(b&&Yd(P)&&(Os(P)?P.value=b[O]:El(P,b[O])),e.state.value[s][O]=P);else if(typeof P=="function"){const U=j(P,O);E[O]=U,p.actions[O]=P}}return ra(x,E),ra(ws(x),E),Object.defineProperty(x,"$state",{get:()=>e.state.value[s],set:O=>{C(P=>{ra(P,O)})}}),e._p.forEach(O=>{ra(x,o.run(()=>O({store:x,app:e._a,pinia:e,options:p})))}),b&&l&&a.hydrate&&a.hydrate(x.$state,b),u=!0,r=!0,x}/*! #__NO_SIDE_EFFECTS__ */function Zr(s,n,a){let e;const t=typeof n=="function";e=t?a:n;function l(o,p){const c=Tr();return o=o||(c?gn(Yr,null):null),o&&Ot(o),o=Xr,o._s.has(s)||(t?Jr(s,n,e,o):Jd(s,e,o)),o._s.get(s)}return l.$id=s,l}const Zd=new Set(["link","style","script","noscript"]),sm=new Set(["title","titleTemplate","script","style","noscript"]),Tl=new Set(["base","meta","link","style","script","noscript"]),nm=new Set(["title","base","htmlAttrs","bodyAttrs","meta","link","style","script","noscript"]),am=new Set(["base","title","titleTemplate","bodyAttrs","htmlAttrs","templateParams"]),em=new Set(["key","tagPosition","tagPriority","tagDuplicateStrategy","innerHTML","textContent","processTemplateParams"]),tm=new Set(["templateParams","htmlAttrs","bodyAttrs"]),lm=new Set(["theme-color","google-site-verification","og","article","book","profile","twitter","author"]);function Al(s,n={},a){for(const e in s){const t=s[e],l=a?`${a}:${e}`:e;typeof t=="object"&&t!==null?Al(t,n,l):typeof t=="function"&&(n[l]=t)}return n}const si=(()=>{if(console.createTask)return console.createTask;const s={run:n=>n()};return()=>s})();function ni(s,n,a,e){for(let t=a;t<s.length;t+=1)try{const l=e?e.run(()=>s[t](...n)):s[t](...n);if(l&&typeof l.then=="function")return Promise.resolve(l).then(()=>ni(s,n,t+1,e))}catch(l){return Promise.reject(l)}}function om(s,n,a){if(s.length>0)return ni(s,n,0,si(a))}function pm(s,n,a){if(s.length>0){const e=si(a);return Promise.all(s.map(t=>e.run(()=>t(...n))))}}function ol(s,n){for(const a of[...s])a(n)}var cm=class{_hooks;_before;_after;_deprecatedHooks;_deprecatedMessages;constructor(){this._hooks={},this._before=void 0,this._after=void 0,this._deprecatedMessages=void 0,this._deprecatedHooks={},this.hook=this.hook.bind(this),this.callHook=this.callHook.bind(this),this.callHookWith=this.callHookWith.bind(this)}hook(s,n,a={}){if(!s||typeof n!="function")return()=>{};const e=s;let t;for(;this._deprecatedHooks[s];)t=this._deprecatedHooks[s],s=t.to;if(t&&!a.allowDeprecated){let l=t.message;l||(l=`${e} hook has been deprecated`+(t.to?`, please use ${t.to}`:"")),this._deprecatedMessages||(this._deprecatedMessages=new Set),this._deprecatedMessages.has(l)||(console.warn(l),this._deprecatedMessages.add(l))}if(!n.name)try{Object.defineProperty(n,"name",{get:()=>"_"+s.replace(/\W+/g,"_")+"_hook_cb",configurable:!0})}catch{}return this._hooks[s]=this._hooks[s]||[],this._hooks[s].push(n),()=>{n&&(this.removeHook(s,n),n=void 0)}}hookOnce(s,n){let a,e=(...t)=>(typeof a=="function"&&a(),a=void 0,e=void 0,n(...t));return a=this.hook(s,e),a}removeHook(s,n){const a=this._hooks[s];if(a){const e=a.indexOf(n);e!==-1&&a.splice(e,1),a.length===0&&(this._hooks[s]=void 0)}}clearHook(s){this._hooks[s]=void 0}deprecateHook(s,n){this._deprecatedHooks[s]=typeof n=="string"?{to:n}:n;const a=this._hooks[s]||[];this._hooks[s]=void 0;for(const e of a)this.hook(s,e)}deprecateHooks(s){for(const n in s)this.deprecateHook(n,s[n])}addHooks(s){const n=Al(s),a=Object.keys(n).map(e=>this.hook(e,n[e]));return()=>{for(const e of a)e();a.length=0}}removeHooks(s){const n=Al(s);for(const a in n)this.removeHook(a,n[a])}removeAllHooks(){this._hooks={}}callHook(s,...n){return this.callHookWith(om,s,n)}callHookParallel(s,...n){return this.callHookWith(pm,s,n)}callHookWith(s,n,a){const e=this._before||this._after?{name:n,args:a,context:{}}:void 0;this._before&&ol(this._before,e);const t=s(this._hooks[n]?[...this._hooks[n]]:[],a,n);return t instanceof Promise?t.finally(()=>{this._after&&e&&ol(this._after,e)}):(this._after&&e&&ol(this._after,e),t)}beforeEach(s){return this._before=this._before||[],this._before.push(s),()=>{if(this._before!==void 0){const n=this._before.indexOf(s);n!==-1&&this._before.splice(n,1)}}}afterEach(s){return this._after=this._after||[],this._after.push(s),()=>{if(this._after!==void 0){const n=this._after.indexOf(s);n!==-1&&this._after.splice(n,1)}}}};function rm(){return new cm}const im=["name","property","http-equiv"],um=new Set(["viewport","description","keywords","robots"]);function ai(s){const n=s.split(":");return n.length?lm.has(n[1]):!1}function Rl(s){const{props:n,tag:a}=s;if(am.has(a))return a;if(a==="link"&&n.rel==="canonical")return"canonical";if(a==="link"&&n.rel==="alternate"){if(n.hreflang)return`alternate:${n.hreflang}`;if(n.type)return`alternate:${n.type}:${n.href||""}`}if(n.charset)return"charset";if(s.tag==="meta"){for(const e of im)if(n[e]!==void 0){const t=n[e],l=t&&typeof t=="string"&&t.includes(":"),o=t&&um.has(t),c=!(l||o)&&s.key?`:key:${s.key}`:"";return`${a}:${t}${c}`}}if(s.key)return`${a}:key:${s.key}`;if(n.id)return`${a}:id:${n.id}`;if(a==="link"&&n.rel==="alternate")return`alternate:${n.href||""}`;if(sm.has(a)){const e=s.textContent||s.innerHTML;if(e)return`${a}:content:${e}`}}function ei(s){const n=s._h||s._d;if(n)return n;const a=s.textContent||s.innerHTML;return a||`${s.tag}:${Object.entries(s.props).map(([e,t])=>`${e}:${String(t)}`).join(",")}`}function vt(s,n,a){typeof s==="function"&&(!a||a!=="titleTemplate"&&!(a[0]==="o"&&a[1]==="n"))&&(s=s());const t=n?n(a,s):s;if(Array.isArray(t))return t.map(l=>vt(l,n));if(t?.constructor===Object){const l={};for(const o of Object.keys(t))l[o]=vt(t[o],n,o);return l}return t}function hm(s,n){const a=s==="style"?new Map:new Set;function e(t){if(t==null||t===void 0)return;const l=String(t).trim();if(l)if(s==="style"){const[o,...p]=l.split(":").map(c=>c?c.trim():"");o&&p.length&&a.set(o,p.join(":"))}else l.split(" ").filter(Boolean).forEach(o=>a.add(o))}return typeof n=="string"?s==="style"?n.split(";").forEach(e):e(n):Array.isArray(n)?n.forEach(t=>e(t)):n&&typeof n=="object"&&Object.entries(n).forEach(([t,l])=>{l&&l!=="false"&&(s==="style"?a.set(String(t).trim(),String(l)):e(t))}),a}function ti(s,n){if(s.props=s.props||{},!n)return s;if(s.tag==="templateParams")return s.props=n,s;const a=Tl.has(s.tag)||s.tag==="htmlAttrs"||s.tag==="bodyAttrs";return Object.entries(n).forEach(([e,t])=>{if(e==="__proto__"||e==="constructor"||e==="prototype")return;if(t===null){s.props[e]=null;return}if(e==="class"||e==="style"){s.props[e]=hm(e,t);return}if(em.has(e)){if((e==="textContent"||e==="innerHTML")&&typeof t=="object"){let u=n.type;if(n.type||(u="application/json"),!u?.endsWith("json")&&u!=="speculationrules")return;n.type=u,s.props.type=u,s[e]=JSON.stringify(t)}else s[e]=t;return}const l=e.startsWith("data-"),o=a&&!l?e.toLowerCase():e,p=String(t),c=s.tag==="meta"&&o==="content";p==="true"||p===""?s.props[o]=l||c?p:!0:!t&&l&&p==="false"?s.props[o]="false":t!==void 0&&(s.props[o]=t)}),s}function dm(s,n){const a=typeof n=="object"&&typeof n!="function"?n:{[s==="script"||s==="noscript"||s==="style"?"innerHTML":"textContent"]:n},e=ti({tag:s,props:{}},a);return e.key&&Zd.has(e.tag)&&(e.props["data-hid"]=e._h=e.key),e.tag==="script"&&typeof e.innerHTML=="object"&&(e.innerHTML=JSON.stringify(e.innerHTML),e.props.type=e.props.type||"application/json"),Array.isArray(e.props.content)?e.props.content.map(t=>({...e,props:{...e.props,content:t}})):e}function mm(s,n){if(!s)return[];typeof s=="function"&&(s=s());const a=(t,l)=>{for(let o=0;o<n.length;o++)l=n[o](t,l);return l};s=a(void 0,s);const e=[];return s=vt(s,a),Object.entries(s||{}).forEach(([t,l])=>{if(l!==void 0)for(const o of Array.isArray(l)?l:[l])e.push(dm(t,o))}),e.flat()}const Ap=(s,n)=>s._w===n._w?s._p-n._p:s._w-n._w,Rp={base:-10,title:10},jm={critical:-8,high:-1,low:2},Lp={meta:{"content-security-policy":-30,charset:-20,viewport:-15},link:{preconnect:20,stylesheet:60,preload:70,modulepreload:70,prefetch:90,"dns-prefetch":90,prerender:90},script:{async:30,defer:80,sync:50},style:{imported:40,sync:60}},gm=/@import/,te=s=>s===""||s===!0;function fm(s,n){if(typeof n.tagPriority=="number")return n.tagPriority;let a=100;const e=jm[n.tagPriority]||0,t=s.resolvedOptions.disableCapoSorting?{link:{},script:{},style:{}}:Lp;if(n.tag in Rp)a=Rp[n.tag];else if(n.tag==="meta"){const l=n.props["http-equiv"]==="content-security-policy"?"content-security-policy":n.props.charset?"charset":n.props.name==="viewport"?"viewport":null;l&&(a=Lp.meta[l])}else if(n.tag==="link"&&n.props.rel)a=t.link[n.props.rel];else if(n.tag==="script"){const l=String(n.props.type);te(n.props.async)?a=t.script.async:n.props.src&&!te(n.props.defer)&&!te(n.props.async)&&l!=="module"&&!l.endsWith("json")||n.innerHTML&&!l.endsWith("json")?a=t.script.sync:(te(n.props.defer)&&n.props.src&&!te(n.props.async)||l==="module")&&(a=t.script.defer)}else n.tag==="style"&&(a=n.innerHTML&&gm.test(n.innerHTML)?t.style.imported:t.style.sync);return(a||100)+e}function Mp(s,n){const a=typeof n=="function"?n(s):n,e=a.key||String(s.plugins.size+1);s.plugins.get(e)||(s.plugins.set(e,a),s.hooks.addHooks(a.hooks||{}))}function bm(s={}){const n=rm();n.addHooks(s.hooks||{});const a=!s.document,e=new Map,t=new Map,l=new Set,o={_entryCount:1,plugins:t,dirty:!1,resolvedOptions:s,hooks:n,ssr:a,entries:e,headEntries(){return[...e.values()]},use:p=>Mp(o,p),push(p,c){const u={...c||{}};delete u.head;const r=u._index??o._entryCount++,i={_i:r,input:p,options:u},h={_poll(d=!1){o.dirty=!0,!d&&l.add(r),n.callHook("entries:updated",o)},dispose(){e.delete(r)&&o.invalidate()},patch(d){(!u.mode||u.mode==="server"&&a||u.mode==="client"&&!a)&&(i.input=d,e.set(r,i),h._poll())}};return h.patch(p),h},async resolveTags(){const p={tagMap:new Map,tags:[],entries:[...o.entries.values()]};for(await n.callHook("entries:resolve",p);l.size;){const h=l.values().next().value;l.delete(h);const d=e.get(h);if(d){const b={tags:mm(d.input,s.propResolvers||[]).map(f=>Object.assign(f,d.options)),entry:d};await n.callHook("entries:normalize",b),d._tags=b.tags.map((f,C)=>(f._w=fm(o,f),f._p=(d._i<<10)+C,f._d=Rl(f),f._d||(f._h=ei(f)),f))}}let c=!1;p.entries.flatMap(h=>(h._tags||[]).map(d=>({...d,props:{...d.props}}))).sort(Ap).reduce((h,d)=>{const b=d._d||d._h;if(!h.has(b))return h.set(b,d);const f=h.get(b);if((d?.tagDuplicateStrategy||(tm.has(d.tag)?"merge":null)||(d.key&&d.key===f.key?"merge":null))==="merge"){const m={...f.props};Object.entries(d.props).forEach(([g,j])=>m[g]=g==="style"?new Map([...f.props.style||new Map,...j]):g==="class"?new Set([...f.props.class||new Set,...j]):j),h.set(b,{...d,props:m})}else d._p>>10===f._p>>10&&d.tag==="meta"&&ai(b)?(h.set(b,Object.assign([...Array.isArray(f)?f:[f],d],d)),c=!0):(d._w===f._w?d._p>f._p:d?._w<f?._w)&&h.set(b,d);return h},p.tagMap);const u=p.tagMap.get("title"),r=p.tagMap.get("titleTemplate");if(o._title=u?.textContent,r){const h=r?.textContent;if(o._titleTemplate=h,h){let d=typeof h=="function"?h(u?.textContent):h;typeof d=="string"&&!o.plugins.has("template-params")&&(d=d.replace("%s",u?.textContent||"")),u?d===null?p.tagMap.delete("title"):p.tagMap.set("title",{...u,textContent:d}):(r.tag="title",r.textContent=d)}}p.tags=Array.from(p.tagMap.values()),c&&(p.tags=p.tags.flat().sort(Ap)),await n.callHook("tags:beforeResolve",p),await n.callHook("tags:resolve",p),await n.callHook("tags:afterResolve",p);const i=[];for(const h of p.tags){const{innerHTML:d,tag:b,props:f}=h;if(nm.has(b)&&!(Object.keys(f).length===0&&!h.innerHTML&&!h.textContent)&&!(b==="meta"&&!f.content&&!f["http-equiv"]&&!f.charset)){if(b==="script"&&d){if(String(f.type).endsWith("json")){const C=typeof d=="string"?d:JSON.stringify(d);h.innerHTML=C.replace(/</g,"\\u003C")}else typeof d=="string"&&(h.innerHTML=d.replace(new RegExp(`</${b}`,"g"),`<\\/${b}`));h._d=Rl(h)}i.push(h)}}return i},invalidate(){for(const p of e.values())l.add(p._i);o.dirty=!0,n.callHook("entries:updated",o)}};return(s?.plugins||[]).forEach(p=>Mp(o,p)),o.hooks.callHook("init",o),s.init?.forEach(p=>p&&o.push(p)),o}const _m=(s,n)=>Os(n)?Yu(n):n,li="usehead";function ym(s){return{install(a){a.config.globalProperties.$unhead=s,a.config.globalProperties.$head=s,a.provide(li,s)}}.install}function vm(){if(Tr()){const s=gn(li);if(s)return s}throw new Error("useHead() was called without provide context, ensure you call it through the setup() function.")}function It(s,n={}){const a=n.head||vm();return a.ssr?a.push(s||{},n):wm(a,s,n)}function wm(s,n,a={}){const e=rs(!1);let t;return Uh(()=>{const o=e.value?{}:vt(n,_m);t?t.patch(o):t=s.push(o,a)}),la()&&(zn(()=>{t.dispose()}),Cr(()=>{e.value=!0}),wr(()=>{e.value=!1})),t}const Cm={created(){let s=!1;const n=la();if(!n)return;const a=n.type;!a||!("head"in a)||(s=typeof a.head=="function"?()=>a.head.call(n.proxy):a.head,s&&It(s))}};async function oi(s,n={}){const a=n.document||s.resolvedOptions.document;if(!a||!s.dirty)return;const e={shouldRender:!0,tags:[]};if(await s.hooks.callHook("dom:beforeRender",e),!!e.shouldRender)return s._domUpdatePromise||(s._domUpdatePromise=new Promise(async t=>{const l=new Map,o=new Promise(d=>{s.resolveTags().then(b=>{d(b.map(f=>{const C=l.get(f._d)||0,m={tag:f,id:(C?`${f._d}:${C}`:f._d)||f._h,shouldRender:!0};return f._d&&ai(f._d)&&l.set(f._d,C+1),m}))})});let p=s._dom;if(!p){p={title:a.title,elMap:new Map().set("htmlAttrs",a.documentElement).set("bodyAttrs",a.body)};for(const d of["body","head"]){const b=a[d]?.children;for(const f of b){const C=f.tagName.toLowerCase();if(!Tl.has(C))continue;const m=ti({tag:C,props:{}},{innerHTML:f.innerHTML,...f.getAttributeNames().reduce((g,j)=>(g[j]=f.getAttribute(j),g),{})||{}});if(m.key=f.getAttribute("data-hid")||void 0,m._d=Rl(m)||ei(m),p.elMap.has(m._d)){let g=1,j=m._d;for(;p.elMap.has(j);)j=`${m._d}:${g++}`;p.elMap.set(j,f)}else p.elMap.set(m._d,f)}}}p.pendingSideEffects={...p.sideEffects},p.sideEffects={};function c(d,b,f){const C=`${d}:${b}`;p.sideEffects[C]=f,delete p.pendingSideEffects[C]}function u({id:d,$el:b,tag:f}){const C=f.tag.endsWith("Attrs");p.elMap.set(d,b),C||(f.textContent&&f.textContent!==b.textContent&&(b.textContent=f.textContent),f.innerHTML&&f.innerHTML!==b.innerHTML&&(b.innerHTML=f.innerHTML),c(d,"el",()=>{b?.remove(),p.elMap.delete(d)}));for(const m in f.props){if(!Object.prototype.hasOwnProperty.call(f.props,m))continue;const g=f.props[m];if(m.startsWith("on")&&typeof g=="function"){const y=b?.dataset;if(y&&y[`${m}fired`]){const x=m.slice(0,-5);g.call(b,new Event(x.substring(2)))}b.getAttribute(`data-${m}`)!==""&&((f.tag==="bodyAttrs"?a.defaultView:b).addEventListener(m.substring(2),g.bind(b)),b.setAttribute(`data-${m}`,""));continue}const j=`attr:${m}`;if(m==="class"){if(!g)continue;for(const y of g)C&&c(d,`${j}:${y}`,()=>b.classList.remove(y)),!b.classList.contains(y)&&b.classList.add(y)}else if(m==="style"){if(!g)continue;for(const[y,x]of g)c(d,`${j}:${y}`,()=>{b.style.removeProperty(y)}),b.style.setProperty(y,x)}else g!==!1&&g!==null&&(b.getAttribute(m)!==g&&b.setAttribute(m,g===!0?"":String(g)),C&&c(d,j,()=>b.removeAttribute(m)))}}const r=[],i={bodyClose:void 0,bodyOpen:void 0,head:void 0},h=await o;for(const d of h){const{tag:b,shouldRender:f,id:C}=d;if(f){if(b.tag==="title"){a.title=b.textContent,c("title","",()=>a.title=p.title);continue}d.$el=d.$el||p.elMap.get(C),d.$el?u(d):Tl.has(b.tag)&&r.push(d)}}for(const d of r){const b=d.tag.tagPosition||"head";d.$el=a.createElement(d.tag.tag),u(d),i[b]=i[b]||a.createDocumentFragment(),i[b].appendChild(d.$el)}for(const d of h)await s.hooks.callHook("dom:renderTag",d,a,c);i.head&&a.head.appendChild(i.head),i.bodyOpen&&a.body.insertBefore(i.bodyOpen,a.body.firstChild),i.bodyClose&&a.body.appendChild(i.bodyClose);for(const d in p.pendingSideEffects)p.pendingSideEffects[d]();s._dom=p,await s.hooks.callHook("dom:rendered",{renders:h}),t()}).finally(()=>{s._domUpdatePromise=void 0,s.dirty=!1})),s._domUpdatePromise}function km(s={}){const n=s.domOptions?.render||oi;s.document=s.document||(typeof window<"u"?document:void 0);const a=s.document?.head.querySelector('script[id="unhead:payload"]')?.innerHTML||!1;return bm({...s,plugins:[...s.plugins||[],{key:"client",hooks:{"entries:updated":n}}],init:[a?JSON.parse(a):!1,...s.init||[]]})}function xm(s,n){let a=0;return()=>{const e=++a;n(()=>{a===e&&s()})}}function Sm(s={}){const n=km({domOptions:{render:xm(()=>oi(n),a=>setTimeout(a,0))},...s});return n.install=ym(n),n}/*!
  * vue-router v4.5.1
  * (c) 2025 Eduardo San Martin Morote
  * @license MIT
  */const za=typeof document<"u";function pi(s){return typeof s=="object"||"displayName"in s||"props"in s||"__vccOpts"in s}function Pm(s){return s.__esModule||s[Symbol.toStringTag]==="Module"||s.default&&pi(s.default)}const ks=Object.assign;function pl(s,n){const a={};for(const e in n){const t=n[e];a[e]=Ln(t)?t.map(s):s(t)}return a}const fe=()=>{},Ln=Array.isArray,ci=/#/g,Em=/&/g,Tm=/\//g,Am=/=/g,Rm=/\?/g,ri=/\+/g,Lm=/%5B/g,Mm=/%5D/g,ii=/%5E/g,Dm=/%60/g,ui=/%7B/g,Om=/%7C/g,hi=/%7D/g,Im=/%20/g;function wo(s){return encodeURI(""+s).replace(Om,"|").replace(Lm,"[").replace(Mm,"]")}function Nm(s){return wo(s).replace(ui,"{").replace(hi,"}").replace(ii,"^")}function Ll(s){return wo(s).replace(ri,"%2B").replace(Im,"+").replace(ci,"%23").replace(Em,"%26").replace(Dm,"`").replace(ui,"{").replace(hi,"}").replace(ii,"^")}function Fm(s){return Ll(s).replace(Am,"%3D")}function $m(s){return wo(s).replace(ci,"%23").replace(Rm,"%3F")}function Bm(s){return s==null?"":$m(s).replace(Tm,"%2F")}function Se(s){try{return decodeURIComponent(""+s)}catch{}return""+s}const qm=/\/$/,zm=s=>s.replace(qm,"");function cl(s,n,a="/"){let e,t={},l="",o="";const p=n.indexOf("#");let c=n.indexOf("?");return p<c&&p>=0&&(c=-1),c>-1&&(e=n.slice(0,c),l=n.slice(c+1,p>-1?p:n.length),t=s(l)),p>-1&&(e=e||n.slice(0,p),o=n.slice(p,n.length)),e=Gm(e??n,a),{fullPath:e+(l&&"?")+l+o,path:e,query:t,hash:Se(o)}}function Vm(s,n){const a=n.query?s(n.query):"";return n.path+(a&&"?")+a+(n.hash||"")}function Dp(s,n){return!n||!s.toLowerCase().startsWith(n.toLowerCase())?s:s.slice(n.length)||"/"}function Um(s,n,a){const e=n.matched.length-1,t=a.matched.length-1;return e>-1&&e===t&&Xa(n.matched[e],a.matched[t])&&di(n.params,a.params)&&s(n.query)===s(a.query)&&n.hash===a.hash}function Xa(s,n){return(s.aliasOf||s)===(n.aliasOf||n)}function di(s,n){if(Object.keys(s).length!==Object.keys(n).length)return!1;for(const a in s)if(!Hm(s[a],n[a]))return!1;return!0}function Hm(s,n){return Ln(s)?Op(s,n):Ln(n)?Op(n,s):s===n}function Op(s,n){return Ln(n)?s.length===n.length&&s.every((a,e)=>a===n[e]):s.length===1&&s[0]===n}function Gm(s,n){if(s.startsWith("/"))return s;if(!s)return n;const a=n.split("/"),e=s.split("/"),t=e[e.length-1];(t===".."||t===".")&&e.push("");let l=a.length-1,o,p;for(o=0;o<e.length;o++)if(p=e[o],p!==".")if(p==="..")l>1&&l--;else break;return a.slice(0,l).join("/")+"/"+e.slice(o).join("/")}const pa={path:"/",name:void 0,params:{},query:{},hash:"",fullPath:"/",matched:[],meta:{},redirectedFrom:void 0};var Pe;(function(s){s.pop="pop",s.push="push"})(Pe||(Pe={}));var be;(function(s){s.back="back",s.forward="forward",s.unknown=""})(be||(be={}));function Wm(s){if(!s)if(za){const n=document.querySelector("base");s=n&&n.getAttribute("href")||"/",s=s.replace(/^\w+:\/\/[^\/]+/,"")}else s="/";return s[0]!=="/"&&s[0]!=="#"&&(s="/"+s),zm(s)}const Km=/^[^#]+#/;function Xm(s,n){return s.replace(Km,"#")+n}function Ym(s,n){const a=document.documentElement.getBoundingClientRect(),e=s.getBoundingClientRect();return{behavior:n.behavior,left:e.left-a.left-(n.left||0),top:e.top-a.top-(n.top||0)}}const Nt=()=>({left:window.scrollX,top:window.scrollY});function Qm(s){let n;if("el"in s){const a=s.el,e=typeof a=="string"&&a.startsWith("#"),t=typeof a=="string"?e?document.getElementById(a.slice(1)):document.querySelector(a):a;if(!t)return;n=Ym(t,s)}else n=s;"scrollBehavior"in document.documentElement.style?window.scrollTo(n):window.scrollTo(n.left!=null?n.left:window.scrollX,n.top!=null?n.top:window.scrollY)}function Ip(s,n){return(history.state?history.state.position-n:-1)+s}const Ml=new Map;function Jm(s,n){Ml.set(s,n)}function Zm(s){const n=Ml.get(s);return Ml.delete(s),n}let sj=()=>location.protocol+"//"+location.host;function mi(s,n){const{pathname:a,search:e,hash:t}=n,l=s.indexOf("#");if(l>-1){let p=t.includes(s.slice(l))?s.slice(l).length:1,c=t.slice(p);return c[0]!=="/"&&(c="/"+c),Dp(c,"")}return Dp(a,s)+e+t}function nj(s,n,a,e){let t=[],l=[],o=null;const p=({state:h})=>{const d=mi(s,location),b=a.value,f=n.value;let C=0;if(h){if(a.value=d,n.value=h,o&&o===b){o=null;return}C=f?h.position-f.position:0}else e(d);t.forEach(m=>{m(a.value,b,{delta:C,type:Pe.pop,direction:C?C>0?be.forward:be.back:be.unknown})})};function c(){o=a.value}function u(h){t.push(h);const d=()=>{const b=t.indexOf(h);b>-1&&t.splice(b,1)};return l.push(d),d}function r(){const{history:h}=window;h.state&&h.replaceState(ks({},h.state,{scroll:Nt()}),"")}function i(){for(const h of l)h();l=[],window.removeEventListener("popstate",p),window.removeEventListener("beforeunload",r)}return window.addEventListener("popstate",p),window.addEventListener("beforeunload",r,{passive:!0}),{pauseListeners:c,listen:u,destroy:i}}function Np(s,n,a,e=!1,t=!1){return{back:s,current:n,forward:a,replaced:e,position:window.history.length,scroll:t?Nt():null}}function aj(s){const{history:n,location:a}=window,e={value:mi(s,a)},t={value:n.state};t.value||l(e.value,{back:null,current:e.value,forward:null,position:n.length-1,replaced:!0,scroll:null},!0);function l(c,u,r){const i=s.indexOf("#"),h=i>-1?(a.host&&document.querySelector("base")?s:s.slice(i))+c:sj()+s+c;try{n[r?"replaceState":"pushState"](u,"",h),t.value=u}catch(d){console.error(d),a[r?"replace":"assign"](h)}}function o(c,u){const r=ks({},n.state,Np(t.value.back,c,t.value.forward,!0),u,{position:t.value.position});l(c,r,!0),e.value=c}function p(c,u){const r=ks({},t.value,n.state,{forward:c,scroll:Nt()});l(r.current,r,!0);const i=ks({},Np(e.value,c,null),{position:r.position+1},u);l(c,i,!1),e.value=c}return{location:e,state:t,push:p,replace:o}}function ej(s){s=Wm(s);const n=aj(s),a=nj(s,n.state,n.location,n.replace);function e(l,o=!0){o||a.pauseListeners(),history.go(l)}const t=ks({location:"",base:s,go:e,createHref:Xm.bind(null,s)},n,a);return Object.defineProperty(t,"location",{enumerable:!0,get:()=>n.location.value}),Object.defineProperty(t,"state",{enumerable:!0,get:()=>n.state.value}),t}function tj(s){return typeof s=="string"||s&&typeof s=="object"}function ji(s){return typeof s=="string"||typeof s=="symbol"}const gi=Symbol("");var Fp;(function(s){s[s.aborted=4]="aborted",s[s.cancelled=8]="cancelled",s[s.duplicated=16]="duplicated"})(Fp||(Fp={}));function Ya(s,n){return ks(new Error,{type:s,[gi]:!0},n)}function Gn(s,n){return s instanceof Error&&gi in s&&(n==null||!!(s.type&n))}const $p="[^/]+?",lj={sensitive:!1,strict:!1,start:!0,end:!0},oj=/[.+*?^${}()[\]/\\]/g;function pj(s,n){const a=ks({},lj,n),e=[];let t=a.start?"^":"";const l=[];for(const u of s){const r=u.length?[]:[90];a.strict&&!u.length&&(t+="/");for(let i=0;i<u.length;i++){const h=u[i];let d=40+(a.sensitive?.25:0);if(h.type===0)i||(t+="/"),t+=h.value.replace(oj,"\\$&"),d+=40;else if(h.type===1){const{value:b,repeatable:f,optional:C,regexp:m}=h;l.push({name:b,repeatable:f,optional:C});const g=m||$p;if(g!==$p){d+=10;try{new RegExp(`(${g})`)}catch(y){throw new Error(`Invalid custom RegExp for param "${b}" (${g}): `+y.message)}}let j=f?`((?:${g})(?:/(?:${g}))*)`:`(${g})`;i||(j=C&&u.length<2?`(?:/${j})`:"/"+j),C&&(j+="?"),t+=j,d+=20,C&&(d+=-8),f&&(d+=-20),g===".*"&&(d+=-50)}r.push(d)}e.push(r)}if(a.strict&&a.end){const u=e.length-1;e[u][e[u].length-1]+=.7000000000000001}a.strict||(t+="/?"),a.end?t+="$":a.strict&&!t.endsWith("/")&&(t+="(?:/|$)");const o=new RegExp(t,a.sensitive?"":"i");function p(u){const r=u.match(o),i={};if(!r)return null;for(let h=1;h<r.length;h++){const d=r[h]||"",b=l[h-1];i[b.name]=d&&b.repeatable?d.split("/"):d}return i}function c(u){let r="",i=!1;for(const h of s){(!i||!r.endsWith("/"))&&(r+="/"),i=!1;for(const d of h)if(d.type===0)r+=d.value;else if(d.type===1){const{value:b,repeatable:f,optional:C}=d,m=b in u?u[b]:"";if(Ln(m)&&!f)throw new Error(`Provided param "${b}" is an array but it is not repeatable (* or + modifiers)`);const g=Ln(m)?m.join("/"):m;if(!g)if(C)h.length<2&&(r.endsWith("/")?r=r.slice(0,-1):i=!0);else throw new Error(`Missing required param "${b}"`);r+=g}}return r||"/"}return{re:o,score:e,keys:l,parse:p,stringify:c}}function cj(s,n){let a=0;for(;a<s.length&&a<n.length;){const e=n[a]-s[a];if(e)return e;a++}return s.length<n.length?s.length===1&&s[0]===80?-1:1:s.length>n.length?n.length===1&&n[0]===80?1:-1:0}function fi(s,n){let a=0;const e=s.score,t=n.score;for(;a<e.length&&a<t.length;){const l=cj(e[a],t[a]);if(l)return l;a++}if(Math.abs(t.length-e.length)===1){if(Bp(e))return 1;if(Bp(t))return-1}return t.length-e.length}function Bp(s){const n=s[s.length-1];return s.length>0&&n[n.length-1]<0}const rj={type:0,value:""},ij=/[a-zA-Z0-9_]/;function uj(s){if(!s)return[[]];if(s==="/")return[[rj]];if(!s.startsWith("/"))throw new Error(`Invalid path "${s}"`);function n(d){throw new Error(`ERR (${a})/"${u}": ${d}`)}let a=0,e=a;const t=[];let l;function o(){l&&t.push(l),l=[]}let p=0,c,u="",r="";function i(){u&&(a===0?l.push({type:0,value:u}):a===1||a===2||a===3?(l.length>1&&(c==="*"||c==="+")&&n(`A repeatable param (${u}) must be alone in its segment. eg: '/:ids+.`),l.push({type:1,value:u,regexp:r,repeatable:c==="*"||c==="+",optional:c==="*"||c==="?"})):n("Invalid state to consume buffer"),u="")}function h(){u+=c}for(;p<s.length;){if(c=s[p++],c==="\\"&&a!==2){e=a,a=4;continue}switch(a){case 0:c==="/"?(u&&i(),o()):c===":"?(i(),a=1):h();break;case 4:h(),a=e;break;case 1:c==="("?a=2:ij.test(c)?h():(i(),a=0,c!=="*"&&c!=="?"&&c!=="+"&&p--);break;case 2:c===")"?r[r.length-1]=="\\"?r=r.slice(0,-1)+c:a=3:r+=c;break;case 3:i(),a=0,c!=="*"&&c!=="?"&&c!=="+"&&p--,r="";break;default:n("Unknown state");break}}return a===2&&n(`Unfinished custom RegExp for param "${u}"`),i(),o(),t}function hj(s,n,a){const e=pj(uj(s.path),a),t=ks(e,{record:s,parent:n,children:[],alias:[]});return n&&!t.record.aliasOf==!n.record.aliasOf&&n.children.push(t),t}function dj(s,n){const a=[],e=new Map;n=Up({strict:!1,end:!0,sensitive:!1},n);function t(i){return e.get(i)}function l(i,h,d){const b=!d,f=zp(i);f.aliasOf=d&&d.record;const C=Up(n,i),m=[f];if("alias"in i){const y=typeof i.alias=="string"?[i.alias]:i.alias;for(const x of y)m.push(zp(ks({},f,{components:d?d.record.components:f.components,path:x,aliasOf:d?d.record:f})))}let g,j;for(const y of m){const{path:x}=y;if(h&&x[0]!=="/"){const R=h.record.path,E=R[R.length-1]==="/"?"":"/";y.path=h.record.path+(x&&E+x)}if(g=hj(y,h,C),d?d.alias.push(g):(j=j||g,j!==g&&j.alias.push(g),b&&i.name&&!Vp(g)&&o(i.name)),bi(g)&&c(g),f.children){const R=f.children;for(let E=0;E<R.length;E++)l(R[E],g,d&&d.children[E])}d=d||g}return j?()=>{o(j)}:fe}function o(i){if(ji(i)){const h=e.get(i);h&&(e.delete(i),a.splice(a.indexOf(h),1),h.children.forEach(o),h.alias.forEach(o))}else{const h=a.indexOf(i);h>-1&&(a.splice(h,1),i.record.name&&e.delete(i.record.name),i.children.forEach(o),i.alias.forEach(o))}}function p(){return a}function c(i){const h=gj(i,a);a.splice(h,0,i),i.record.name&&!Vp(i)&&e.set(i.record.name,i)}function u(i,h){let d,b={},f,C;if("name"in i&&i.name){if(d=e.get(i.name),!d)throw Ya(1,{location:i});C=d.record.name,b=ks(qp(h.params,d.keys.filter(j=>!j.optional).concat(d.parent?d.parent.keys.filter(j=>j.optional):[]).map(j=>j.name)),i.params&&qp(i.params,d.keys.map(j=>j.name))),f=d.stringify(b)}else if(i.path!=null)f=i.path,d=a.find(j=>j.re.test(f)),d&&(b=d.parse(f),C=d.record.name);else{if(d=h.name?e.get(h.name):a.find(j=>j.re.test(h.path)),!d)throw Ya(1,{location:i,currentLocation:h});C=d.record.name,b=ks({},h.params,i.params),f=d.stringify(b)}const m=[];let g=d;for(;g;)m.unshift(g.record),g=g.parent;return{name:C,path:f,params:b,matched:m,meta:jj(m)}}s.forEach(i=>l(i));function r(){a.length=0,e.clear()}return{addRoute:l,resolve:u,removeRoute:o,clearRoutes:r,getRoutes:p,getRecordMatcher:t}}function qp(s,n){const a={};for(const e of n)e in s&&(a[e]=s[e]);return a}function zp(s){const n={path:s.path,redirect:s.redirect,name:s.name,meta:s.meta||{},aliasOf:s.aliasOf,beforeEnter:s.beforeEnter,props:mj(s),children:s.children||[],instances:{},leaveGuards:new Set,updateGuards:new Set,enterCallbacks:{},components:"components"in s?s.components||null:s.component&&{default:s.component}};return Object.defineProperty(n,"mods",{value:{}}),n}function mj(s){const n={},a=s.props||!1;if("component"in s)n.default=a;else for(const e in s.components)n[e]=typeof a=="object"?a[e]:a;return n}function Vp(s){for(;s;){if(s.record.aliasOf)return!0;s=s.parent}return!1}function jj(s){return s.reduce((n,a)=>ks(n,a.meta),{})}function Up(s,n){const a={};for(const e in s)a[e]=e in n?n[e]:s[e];return a}function gj(s,n){let a=0,e=n.length;for(;a!==e;){const l=a+e>>1;fi(s,n[l])<0?e=l:a=l+1}const t=fj(s);return t&&(e=n.lastIndexOf(t,e-1)),e}function fj(s){let n=s;for(;n=n.parent;)if(bi(n)&&fi(s,n)===0)return n}function bi({record:s}){return!!(s.name||s.components&&Object.keys(s.components).length||s.redirect)}function bj(s){const n={};if(s===""||s==="?")return n;const e=(s[0]==="?"?s.slice(1):s).split("&");for(let t=0;t<e.length;++t){const l=e[t].replace(ri," "),o=l.indexOf("="),p=Se(o<0?l:l.slice(0,o)),c=o<0?null:Se(l.slice(o+1));if(p in n){let u=n[p];Ln(u)||(u=n[p]=[u]),u.push(c)}else n[p]=c}return n}function Hp(s){let n="";for(let a in s){const e=s[a];if(a=Fm(a),e==null){e!==void 0&&(n+=(n.length?"&":"")+a);continue}(Ln(e)?e.map(l=>l&&Ll(l)):[e&&Ll(e)]).forEach(l=>{l!==void 0&&(n+=(n.length?"&":"")+a,l!=null&&(n+="="+l))})}return n}function _j(s){const n={};for(const a in s){const e=s[a];e!==void 0&&(n[a]=Ln(e)?e.map(t=>t==null?null:""+t):e==null?e:""+e)}return n}const yj=Symbol(""),Gp=Symbol(""),Ft=Symbol(""),Co=Symbol(""),Dl=Symbol("");function le(){let s=[];function n(e){return s.push(e),()=>{const t=s.indexOf(e);t>-1&&s.splice(t,1)}}function a(){s=[]}return{add:n,list:()=>s.slice(),reset:a}}function ha(s,n,a,e,t,l=o=>o()){const o=e&&(e.enterCallbacks[t]=e.enterCallbacks[t]||[]);return()=>new Promise((p,c)=>{const u=h=>{h===!1?c(Ya(4,{from:a,to:n})):h instanceof Error?c(h):tj(h)?c(Ya(2,{from:n,to:h})):(o&&e.enterCallbacks[t]===o&&typeof h=="function"&&o.push(h),p())},r=l(()=>s.call(e&&e.instances[t],n,a,u));let i=Promise.resolve(r);s.length<3&&(i=i.then(u)),i.catch(h=>c(h))})}function rl(s,n,a,e,t=l=>l()){const l=[];for(const o of s)for(const p in o.components){let c=o.components[p];if(!(n!=="beforeRouteEnter"&&!o.instances[p]))if(pi(c)){const r=(c.__vccOpts||c)[n];r&&l.push(ha(r,a,e,o,p,t))}else{let u=c();l.push(()=>u.then(r=>{if(!r)throw new Error(`Couldn't resolve component "${p}" at "${o.path}"`);const i=Pm(r)?r.default:r;o.mods[p]=r,o.components[p]=i;const d=(i.__vccOpts||i)[n];return d&&ha(d,a,e,o,p,t)()}))}}return l}function Wp(s){const n=gn(Ft),a=gn(Co),e=ns(()=>{const c=ss(s.to);return n.resolve(c)}),t=ns(()=>{const{matched:c}=e.value,{length:u}=c,r=c[u-1],i=a.matched;if(!r||!i.length)return-1;const h=i.findIndex(Xa.bind(null,r));if(h>-1)return h;const d=Kp(c[u-2]);return u>1&&Kp(r)===d&&i[i.length-1].path!==d?i.findIndex(Xa.bind(null,c[u-2])):h}),l=ns(()=>t.value>-1&&xj(a.params,e.value.params)),o=ns(()=>t.value>-1&&t.value===a.matched.length-1&&di(a.params,e.value.params));function p(c={}){if(kj(c)){const u=n[ss(s.replace)?"replace":"push"](ss(s.to)).catch(fe);return s.viewTransition&&typeof document<"u"&&"startViewTransition"in document&&document.startViewTransition(()=>u),u}return Promise.resolve()}return{route:e,href:ns(()=>e.value.href),isActive:l,isExactActive:o,navigate:p}}function vj(s){return s.length===1?s[0]:s}const wj=Ns({name:"RouterLink",compatConfig:{MODE:3},props:{to:{type:[String,Object],required:!0},replace:Boolean,activeClass:String,exactActiveClass:String,custom:Boolean,ariaCurrentValue:{type:String,default:"page"},viewTransition:Boolean},useLink:Wp,setup(s,{slots:n}){const a=Me(Wp(s)),{options:e}=gn(Ft),t=ns(()=>({[Xp(s.activeClass,e.linkActiveClass,"router-link-active")]:a.isActive,[Xp(s.exactActiveClass,e.linkExactActiveClass,"router-link-exact-active")]:a.isExactActive}));return()=>{const l=n.default&&vj(n.default(a));return s.custom?l:Be("a",{"aria-current":a.isExactActive?s.ariaCurrentValue:null,href:a.href,onClick:a.navigate,class:t.value},l)}}}),Cj=wj;function kj(s){if(!(s.metaKey||s.altKey||s.ctrlKey||s.shiftKey)&&!s.defaultPrevented&&!(s.button!==void 0&&s.button!==0)){if(s.currentTarget&&s.currentTarget.getAttribute){const n=s.currentTarget.getAttribute("target");if(/\b_blank\b/i.test(n))return}return s.preventDefault&&s.preventDefault(),!0}}function xj(s,n){for(const a in n){const e=n[a],t=s[a];if(typeof e=="string"){if(e!==t)return!1}else if(!Ln(t)||t.length!==e.length||e.some((l,o)=>l!==t[o]))return!1}return!0}function Kp(s){return s?s.aliasOf?s.aliasOf.path:s.path:""}const Xp=(s,n,a)=>s??n??a,Sj=Ns({name:"RouterView",inheritAttrs:!1,props:{name:{type:String,default:"default"},route:Object},compatConfig:{MODE:3},setup(s,{attrs:n,slots:a}){const e=gn(Dl),t=ns(()=>s.route||e.value),l=gn(Gp,0),o=ns(()=>{let u=ss(l);const{matched:r}=t.value;let i;for(;(i=r[u])&&!i.components;)u++;return u}),p=ns(()=>t.value.matched[o.value]);et(Gp,ns(()=>o.value+1)),et(yj,p),et(Dl,t);const c=rs();return Vs(()=>[c.value,p.value,s.name],([u,r,i],[h,d,b])=>{r&&(r.instances[i]=u,d&&d!==r&&u&&u===h&&(r.leaveGuards.size||(r.leaveGuards=d.leaveGuards),r.updateGuards.size||(r.updateGuards=d.updateGuards))),u&&r&&(!d||!Xa(r,d)||!h)&&(r.enterCallbacks[i]||[]).forEach(f=>f(u))},{flush:"post"}),()=>{const u=t.value,r=s.name,i=p.value,h=i&&i.components[r];if(!h)return Yp(a.default,{Component:h,route:u});const d=i.props[r],b=d?d===!0?u.params:typeof d=="function"?d(u):d:null,C=Be(h,ks({},b,n,{onVnodeUnmounted:m=>{m.component.isUnmounted&&(i.instances[r]=null)},ref:c}));return Yp(a.default,{Component:C,route:u})||C}}});function Yp(s,n){if(!s)return null;const a=s(n);return a.length===1?a[0]:a}const Pj=Sj;function Ej(s){const n=dj(s.routes,s),a=s.parseQuery||bj,e=s.stringifyQuery||Hp,t=s.history,l=le(),o=le(),p=le(),c=io(pa);let u=pa;za&&s.scrollBehavior&&"scrollRestoration"in history&&(history.scrollRestoration="manual");const r=pl.bind(null,L=>""+L),i=pl.bind(null,Bm),h=pl.bind(null,Se);function d(L,W){let H,J;return ji(L)?(H=n.getRecordMatcher(L),J=W):J=L,n.addRoute(J,H)}function b(L){const W=n.getRecordMatcher(L);W&&n.removeRoute(W)}function f(){return n.getRoutes().map(L=>L.record)}function C(L){return!!n.getRecordMatcher(L)}function m(L,W){if(W=ks({},W||c.value),typeof L=="string"){const A=cl(a,L,W.path),F=n.resolve({path:A.path},W),q=t.createHref(A.fullPath);return ks(A,F,{params:h(F.params),hash:Se(A.hash),redirectedFrom:void 0,href:q})}let H;if(L.path!=null)H=ks({},L,{path:cl(a,L.path,W.path).path});else{const A=ks({},L.params);for(const F in A)A[F]==null&&delete A[F];H=ks({},L,{params:i(A)}),W.params=i(W.params)}const J=n.resolve(H,W),hs=L.hash||"";J.params=r(h(J.params));const w=Vm(e,ks({},L,{hash:Nm(hs),path:J.path})),k=t.createHref(w);return ks({fullPath:w,hash:hs,query:e===Hp?_j(L.query):L.query||{}},J,{redirectedFrom:void 0,href:k})}function g(L){return typeof L=="string"?cl(a,L,c.value.path):ks({},L)}function j(L,W){if(u!==L)return Ya(8,{from:W,to:L})}function y(L){return E(L)}function x(L){return y(ks(g(L),{replace:!0}))}function R(L){const W=L.matched[L.matched.length-1];if(W&&W.redirect){const{redirect:H}=W;let J=typeof H=="function"?H(L):H;return typeof J=="string"&&(J=J.includes("?")||J.includes("#")?J=g(J):{path:J},J.params={}),ks({query:L.query,hash:L.hash,params:J.path!=null?{}:L.params},J)}}function E(L,W){const H=u=m(L),J=c.value,hs=L.state,w=L.force,k=L.replace===!0,A=R(H);if(A)return E(ks(g(A),{state:typeof A=="object"?ks({},hs,A.state):hs,force:w,replace:k}),W||H);const F=H;F.redirectedFrom=W;let q;return!w&&Um(e,J,H)&&(q=Ya(16,{to:F,from:J}),Ms(J,J,!0,!1)),(q?Promise.resolve(q):U(F,J)).catch(B=>Gn(B)?Gn(B,2)?B:Qs(B):us(B,F,J)).then(B=>{if(B){if(Gn(B,2))return E(ks({replace:k},g(B.to),{state:typeof B.to=="object"?ks({},hs,B.to.state):hs,force:w}),W||F)}else B=V(F,J,!0,k,hs);return Y(F,J,B),B})}function O(L,W){const H=j(L,W);return H?Promise.reject(H):Promise.resolve()}function P(L){const W=en.values().next().value;return W&&typeof W.runWithContext=="function"?W.runWithContext(L):L()}function U(L,W){let H;const[J,hs,w]=Tj(L,W);H=rl(J.reverse(),"beforeRouteLeave",L,W);for(const A of J)A.leaveGuards.forEach(F=>{H.push(ha(F,L,W))});const k=O.bind(null,L,W);return H.push(k),os(H).then(()=>{H=[];for(const A of l.list())H.push(ha(A,L,W));return H.push(k),os(H)}).then(()=>{H=rl(hs,"beforeRouteUpdate",L,W);for(const A of hs)A.updateGuards.forEach(F=>{H.push(ha(F,L,W))});return H.push(k),os(H)}).then(()=>{H=[];for(const A of w)if(A.beforeEnter)if(Ln(A.beforeEnter))for(const F of A.beforeEnter)H.push(ha(F,L,W));else H.push(ha(A.beforeEnter,L,W));return H.push(k),os(H)}).then(()=>(L.matched.forEach(A=>A.enterCallbacks={}),H=rl(w,"beforeRouteEnter",L,W,P),H.push(k),os(H))).then(()=>{H=[];for(const A of o.list())H.push(ha(A,L,W));return H.push(k),os(H)}).catch(A=>Gn(A,8)?A:Promise.reject(A))}function Y(L,W,H){p.list().forEach(J=>P(()=>J(L,W,H)))}function V(L,W,H,J,hs){const w=j(L,W);if(w)return w;const k=W===pa,A=za?history.state:{};H&&(J||k?t.replace(L.fullPath,ks({scroll:k&&A&&A.scroll},hs)):t.push(L.fullPath,hs)),c.value=L,Ms(L,W,H,k),Qs()}let z;function Q(){z||(z=t.listen((L,W,H)=>{if(!Mn.listening)return;const J=m(L),hs=R(J);if(hs){E(ks(hs,{replace:!0,force:!0}),J).catch(fe);return}u=J;const w=c.value;za&&Jm(Ip(w.fullPath,H.delta),Nt()),U(J,w).catch(k=>Gn(k,12)?k:Gn(k,2)?(E(ks(g(k.to),{force:!0}),J).then(A=>{Gn(A,20)&&!H.delta&&H.type===Pe.pop&&t.go(-1,!1)}).catch(fe),Promise.reject()):(H.delta&&t.go(-H.delta,!1),us(k,J,w))).then(k=>{k=k||V(J,w,!1),k&&(H.delta&&!Gn(k,8)?t.go(-H.delta,!1):H.type===Pe.pop&&Gn(k,20)&&t.go(-1,!1)),Y(J,w,k)}).catch(fe)}))}let cs=le(),ls=le(),X;function us(L,W,H){Qs(L);const J=ls.list();return J.length?J.forEach(hs=>hs(L,W,H)):console.error(L),Promise.reject(L)}function Fs(){return X&&c.value!==pa?Promise.resolve():new Promise((L,W)=>{cs.add([L,W])})}function Qs(L){return X||(X=!L,Q(),cs.list().forEach(([W,H])=>L?H(L):W()),cs.reset()),L}function Ms(L,W,H,J){const{scrollBehavior:hs}=s;if(!za||!hs)return Promise.resolve();const w=!H&&Zm(Ip(L.fullPath,0))||(J||!H)&&history.state&&history.state.scroll||null;return bn().then(()=>hs(L,W,w)).then(k=>k&&Qm(k)).catch(k=>us(k,L,W))}const $s=L=>t.go(L);let dn;const en=new Set,Mn={currentRoute:c,listening:!0,addRoute:d,removeRoute:b,clearRoutes:n.clearRoutes,hasRoute:C,getRoutes:f,resolve:m,options:s,push:y,replace:x,go:$s,back:()=>$s(-1),forward:()=>$s(1),beforeEach:l.add,beforeResolve:o.add,afterEach:p.add,onError:ls.add,isReady:Fs,install(L){const W=this;L.component("RouterLink",Cj),L.component("RouterView",Pj),L.config.globalProperties.$router=W,Object.defineProperty(L.config.globalProperties,"$route",{enumerable:!0,get:()=>ss(c)}),za&&!dn&&c.value===pa&&(dn=!0,y(t.location).catch(hs=>{}));const H={};for(const hs in pa)Object.defineProperty(H,hs,{get:()=>c.value[hs],enumerable:!0});L.provide(Ft,W),L.provide(Co,lr(H)),L.provide(Dl,c);const J=L.unmount;en.add(L),L.unmount=function(){en.delete(L),en.size<1&&(u=pa,z&&z(),z=null,c.value=pa,dn=!1,X=!1),J()}}};function os(L){return L.reduce((W,H)=>W.then(()=>P(H)),Promise.resolve())}return Mn}function Tj(s,n){const a=[],e=[],t=[],l=Math.max(n.matched.length,s.matched.length);for(let o=0;o<l;o++){const p=n.matched[o];p&&(s.matched.find(u=>Xa(u,p))?e.push(p):a.push(p));const c=s.matched[o];c&&(n.matched.find(u=>Xa(u,c))||t.push(c))}return[a,e,t]}function qe(){return gn(Ft)}function ya(s){return gn(Co)}function Aj(s){return document.readyState==="loading"?new Promise(n=>{document.addEventListener("DOMContentLoaded",()=>n(s))}):Promise.resolve(s)}const Rj=Ns({setup(s,{slots:n}){const a=rs(!1);return En(()=>a.value=!0),()=>a.value?n.default&&n.default({}):n.placeholder&&n.placeholder({})}});function Lj(s){try{return JSON.parse(s||"{}")}catch(n){return console.error("[SSG] On state deserialization -",n,s),{}}}function Mj(s,n,a,e){const{transformState:t,registerComponents:l=!0,useHead:o=!0,rootContainer:p="#app"}={};async function c(u){const r=Ud(s);let i;o&&r.use(i=Sm());const h=Ej({history:ej(n.base),...n}),{routes:d}=n;l&&r.component("ClientOnly",Rj);const b=[],m={app:r,head:i,isClient:!0,router:h,routes:d,onSSRAppRendered:()=>{},triggerOnSSRAppRendered:()=>Promise.all(b.map(x=>x())),initialState:{},transformState:t,routePath:u};await Aj(),m.initialState=t?.(window.__INITIAL_STATE__||{})||Lj(window.__INITIAL_STATE__),await a?.(m),r.use(h);let g,j=!0;h.beforeEach((x,R,E)=>{(j||g&&g===x.path)&&(j=!1,g=x.path,x.meta.state=m.initialState),E()});const y=m.initialState;return{...m,initialState:y}}return(async()=>{const{app:u,router:r}=await c();await r.isReady(),u.mount(p,!0)})(),c}/*!
  * vue-i18n v12.0.0-alpha.3
  * (c) 2016-present kazuya kawaguchi and contributors
  * Released under the MIT License.
  */const aa=typeof window<"u";let wn,Da;{const s=aa&&window.performance;s&&s.mark&&s.measure&&s.clearMarks&&s.clearMeasures&&(wn=n=>{s.mark(n)},Da=(n,a,e)=>{s.measure(n,a,e),s.clearMarks(a),s.clearMarks(e)})}const Dj=/\{([0-9a-zA-Z]+)\}/g;function $t(s,...n){return n.length===1&&vs(n[0])&&(n=n[0]),(!n||!n.hasOwnProperty)&&(n={}),s.replace(Dj,(a,e)=>n.hasOwnProperty(e)?n[e]:"")}const Vn=(s,n=!1)=>n?Symbol.for(s):Symbol(s),Oj=(s,n,a)=>Ij({l:s,k:n,s:a}),Ij=s=>JSON.stringify(s).replace(/\u2028/g,"\\u2028").replace(/\u2029/g,"\\u2029").replace(/\u0027/g,"\\u0027"),Hs=s=>typeof s=="number"&&isFinite(s),Nj=s=>ko(s)==="[object Date]",wt=s=>ko(s)==="[object RegExp]",Bt=s=>ys(s)&&Object.keys(s).length===0,Ys=Object.assign,Fj=Object.create,Ps=(s=null)=>Fj(s);let Qp;const $j=()=>Qp||(Qp=typeof globalThis<"u"?globalThis:typeof self<"u"?self:typeof window<"u"?window:typeof global<"u"?global:Ps());function Jp(s){return s.replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;").replace(/'/g,"&apos;")}const Bj=Object.prototype.hasOwnProperty;function Tn(s,n){return Bj.call(s,n)}const qs=Array.isArray,Ls=s=>typeof s=="function",ts=s=>typeof s=="string",Ds=s=>typeof s=="boolean",vs=s=>s!==null&&typeof s=="object",qj=s=>vs(s)&&Ls(s.then)&&Ls(s.catch),_i=Object.prototype.toString,ko=s=>_i.call(s),ys=s=>ko(s)==="[object Object]",zj=s=>s==null?"":qs(s)||ys(s)&&s.toString===_i?JSON.stringify(s,null,2):String(s);function yi(s,n=""){return s.reduce((a,e,t)=>t===0?a+e:a+n+e,"")}const Zp=2;function Vj(s,n=0,a=s.length){const e=s.split(/\r?\n/);let t=0;const l=[];for(let o=0;o<e.length;o++)if(t+=e[o].length+1,t>=n){for(let p=o-Zp;p<=o+Zp||a>t;p++){if(p<0||p>=e.length)continue;const c=p+1;l.push(`${c}${" ".repeat(3-String(c).length)}|  ${e[p]}`);const u=e[p].length;if(p===o){const r=n-(t-u)+1,i=Math.max(1,a>t?u-r:a-n);l.push("   |  "+" ".repeat(r)+"^".repeat(i))}else if(p>o){if(a>t){const r=Math.max(Math.min(a-t,u),1);l.push("   |  "+"^".repeat(r))}t+=u+1}}break}return l.join(`
`)}function va(s,n){typeof console<"u"&&(console.warn("[intlify] "+s),n&&console.warn(n.stack))}const sc={};function Uj(s){sc[s]||(sc[s]=!0,va(s))}function vi(){const s=new Map;return{events:s,on(a,e){const t=s.get(a);t&&t.push(e)||s.set(a,[e])},off(a,e){const t=s.get(a);t&&t.splice(t.indexOf(e)>>>0,1)},emit(a,e){(s.get(a)||[]).slice().map(t=>t(e)),(s.get("*")||[]).slice().map(t=>t(a,e))}}}const Ye=s=>!vs(s)||qs(s);function pt(s,n){if(Ye(s)||Ye(n))throw new Error("Invalid value");const a=[{src:s,des:n}];for(;a.length;){const{src:e,des:t}=a.pop();Object.keys(e).forEach(l=>{l!=="__proto__"&&(vs(e[l])&&!vs(t[l])&&(t[l]=Array.isArray(e[l])?[]:Ps()),Ye(t[l])||Ye(e[l])?t[l]=e[l]:a.push({src:e[l],des:t[l]}))})}}function Hj(s,n,a){return{line:s,column:n,offset:a}}function Ol(s,n,a){return{start:s,end:n}}const ds={EXPECTED_TOKEN:1,INVALID_TOKEN_IN_PLACEHOLDER:2,UNTERMINATED_SINGLE_QUOTE_IN_PLACEHOLDER:3,UNKNOWN_ESCAPE_SEQUENCE:4,INVALID_UNICODE_ESCAPE_SEQUENCE:5,UNBALANCED_CLOSING_BRACE:6,UNTERMINATED_CLOSING_BRACE:7,EMPTY_PLACEHOLDER:8,NOT_ALLOW_NEST_PLACEHOLDER:9,INVALID_LINKED_FORMAT:10,MUST_HAVE_MESSAGES_IN_PLURAL:11,UNEXPECTED_EMPTY_LINKED_MODIFIER:12,UNEXPECTED_EMPTY_LINKED_KEY:13,UNEXPECTED_LEXICAL_ANALYSIS:14,UNHANDLED_CODEGEN_NODE_TYPE:15,UNHANDLED_MINIFIER_NODE_TYPE:16},Gj=17,Wj={[ds.EXPECTED_TOKEN]:"Expected token: '{0}'",[ds.INVALID_TOKEN_IN_PLACEHOLDER]:"Invalid token in placeholder: '{0}'",[ds.UNTERMINATED_SINGLE_QUOTE_IN_PLACEHOLDER]:"Unterminated single quote in placeholder",[ds.UNKNOWN_ESCAPE_SEQUENCE]:"Unknown escape sequence: \\{0}",[ds.INVALID_UNICODE_ESCAPE_SEQUENCE]:"Invalid unicode escape sequence: {0}",[ds.UNBALANCED_CLOSING_BRACE]:"Unbalanced closing brace",[ds.UNTERMINATED_CLOSING_BRACE]:"Unterminated closing brace",[ds.EMPTY_PLACEHOLDER]:"Empty placeholder",[ds.NOT_ALLOW_NEST_PLACEHOLDER]:"Not allowed nest placeholder",[ds.INVALID_LINKED_FORMAT]:"Invalid linked format",[ds.MUST_HAVE_MESSAGES_IN_PLURAL]:"Plural must have messages",[ds.UNEXPECTED_EMPTY_LINKED_MODIFIER]:"Unexpected empty linked modifier",[ds.UNEXPECTED_EMPTY_LINKED_KEY]:"Unexpected empty linked key",[ds.UNEXPECTED_LEXICAL_ANALYSIS]:"Unexpected lexical analysis in token: '{0}'",[ds.UNHANDLED_CODEGEN_NODE_TYPE]:"unhandled codegen node type: '{0}'",[ds.UNHANDLED_MINIFIER_NODE_TYPE]:"unhandled mimifier node type: '{0}'"};function ze(s,n,a={}){const{domain:e,messages:t,args:l}=a,o=$t((t||Wj)[s]||"",...l||[]),p=new SyntaxError(String(o));return p.code=s,n&&(p.location=n),p.domain=e,p}function Kj(s){throw s}const Xj="minifier";function Va(s){switch(s.t=s.type,s.type){case 0:{const n=s;Va(n.body),n.b=n.body,delete n.body;break}case 1:{const n=s,a=n.cases;for(let e=0;e<a.length;e++)Va(a[e]);n.c=a,delete n.cases;break}case 2:{const n=s,a=n.items;for(let e=0;e<a.length;e++)Va(a[e]);n.i=a,delete n.items,n.static&&(n.s=n.static,delete n.static);break}case 3:case 9:case 8:case 7:{const n=s;n.value&&(n.v=n.value,delete n.value);break}case 6:{const n=s;Va(n.key),n.k=n.key,delete n.key,n.modifier&&(Va(n.modifier),n.m=n.modifier,delete n.modifier);break}case 5:{const n=s;n.i=n.index,delete n.index;break}case 4:{const n=s;n.k=n.key,delete n.key;break}default:throw ze(ds.UNHANDLED_MINIFIER_NODE_TYPE,null,{domain:Xj,args:[s.type]})}delete s.type}function Yj(s){const n=s.body;return n.type===2?nc(n):n.cases.forEach(a=>nc(a)),s}function nc(s){if(s.items.length===1){const n=s.items[0];(n.type===3||n.type===9)&&(s.static=n.value,delete n.value)}else{const n=[];for(let a=0;a<s.items.length;a++){const e=s.items[a];if(!(e.type===3||e.type===9)||e.value==null)break;n.push(e.value)}if(n.length===s.items.length){s.static=yi(n);for(let a=0;a<s.items.length;a++){const e=s.items[a];(e.type===3||e.type===9)&&delete e.value}}}}const Wn=" ",Qj="\r",on=`
`,Jj="\u2028",Zj="\u2029";function sg(s){const n=s;let a=0,e=1,t=1,l=0;const o=E=>n[E]===Qj&&n[E+1]===on,p=E=>n[E]===on,c=E=>n[E]===Zj,u=E=>n[E]===Jj,r=E=>o(E)||p(E)||c(E)||u(E),i=()=>a,h=()=>e,d=()=>t,b=()=>l,f=E=>o(E)||c(E)||u(E)?on:n[E],C=()=>f(a),m=()=>f(a+l);function g(){return l=0,r(a)&&(e++,t=0),o(a)&&a++,a++,t++,n[a]}function j(){return o(a+l)&&l++,l++,n[a+l]}function y(){a=0,e=1,t=1,l=0}function x(E=0){l=E}function R(){const E=a+l;for(;E!==a;)g();l=0}return{index:i,line:h,column:d,peekOffset:b,charAt:f,currentChar:C,currentPeek:m,next:g,peek:j,reset:y,resetPeek:x,skipToPeek:R}}const ca=void 0,ng=".",ac="'",ag="tokenizer";function eg(s,n={}){const a=n.location!==!1,e=sg(s),t=()=>e.index(),l=()=>Hj(e.line(),e.column(),e.index()),o=l(),p=t(),c={currentType:13,offset:p,startLoc:o,endLoc:o,lastType:13,lastOffset:p,lastStartLoc:o,lastEndLoc:o,braceNest:0,inLinked:!1,text:""},u=()=>c,{onError:r}=n;function i(_,v,T,...D){const Z=u();if(v.column+=T,v.offset+=T,r){const K=a?Ol(Z.startLoc,v):null,as=ze(_,K,{domain:ag,args:D});r(as)}}function h(_,v,T){_.endLoc=l(),_.currentType=v;const D={type:v};return a&&(D.loc=Ol(_.startLoc,_.endLoc)),T!=null&&(D.value=T),D}const d=_=>h(_,13);function b(_,v){return _.currentChar()===v?(_.next(),v):(i(ds.EXPECTED_TOKEN,l(),0,v),"")}function f(_){let v="";for(;_.currentPeek()===Wn||_.currentPeek()===on;)v+=_.currentPeek(),_.peek();return v}function C(_){const v=f(_);return _.skipToPeek(),v}function m(_){if(_===ca)return!1;const v=_.charCodeAt(0);return v>=97&&v<=122||v>=65&&v<=90||v===95}function g(_){if(_===ca)return!1;const v=_.charCodeAt(0);return v>=48&&v<=57}function j(_,v){const{currentType:T}=v;if(T!==2)return!1;f(_);const D=m(_.currentPeek());return _.resetPeek(),D}function y(_,v){const{currentType:T}=v;if(T!==2)return!1;f(_);const D=_.currentPeek()==="-"?_.peek():_.currentPeek(),Z=g(D);return _.resetPeek(),Z}function x(_,v){const{currentType:T}=v;if(T!==2)return!1;f(_);const D=_.currentPeek()===ac;return _.resetPeek(),D}function R(_,v){const{currentType:T}=v;if(T!==7)return!1;f(_);const D=_.currentPeek()===".";return _.resetPeek(),D}function E(_,v){const{currentType:T}=v;if(T!==8)return!1;f(_);const D=m(_.currentPeek());return _.resetPeek(),D}function O(_,v){const{currentType:T}=v;if(!(T===7||T===11))return!1;f(_);const D=_.currentPeek()===":";return _.resetPeek(),D}function P(_,v){const{currentType:T}=v;if(T!==9)return!1;const D=()=>{const K=_.currentPeek();return K==="{"?m(_.peek()):K==="@"||K==="|"||K===":"||K==="."||K===Wn||!K?!1:K===on?(_.peek(),D()):Y(_,!1)},Z=D();return _.resetPeek(),Z}function U(_){f(_);const v=_.currentPeek()==="|";return _.resetPeek(),v}function Y(_,v=!0){const T=(Z=!1,K="")=>{const as=_.currentPeek();return as==="{"||as==="@"||!as?Z:as==="|"?!(K===Wn||K===on):as===Wn?(_.peek(),T(!0,Wn)):as===on?(_.peek(),T(!0,on)):!0},D=T();return v&&_.resetPeek(),D}function V(_,v){const T=_.currentChar();return T===ca?ca:v(T)?(_.next(),T):null}function z(_){const v=_.charCodeAt(0);return v>=97&&v<=122||v>=65&&v<=90||v>=48&&v<=57||v===95||v===36}function Q(_){return V(_,z)}function cs(_){const v=_.charCodeAt(0);return v>=97&&v<=122||v>=65&&v<=90||v>=48&&v<=57||v===95||v===36||v===45}function ls(_){return V(_,cs)}function X(_){const v=_.charCodeAt(0);return v>=48&&v<=57}function us(_){return V(_,X)}function Fs(_){const v=_.charCodeAt(0);return v>=48&&v<=57||v>=65&&v<=70||v>=97&&v<=102}function Qs(_){return V(_,Fs)}function Ms(_){let v="",T="";for(;v=us(_);)T+=v;return T}function $s(_){let v="";for(;;){const T=_.currentChar();if(T==="{"||T==="}"||T==="@"||T==="|"||!T)break;if(T===Wn||T===on)if(Y(_))v+=T,_.next();else{if(U(_))break;v+=T,_.next()}else v+=T,_.next()}return v}function dn(_){C(_);let v="",T="";for(;v=ls(_);)T+=v;return _.currentChar()===ca&&i(ds.UNTERMINATED_CLOSING_BRACE,l(),0),T}function en(_){C(_);let v="";return _.currentChar()==="-"?(_.next(),v+=`-${Ms(_)}`):v+=Ms(_),_.currentChar()===ca&&i(ds.UNTERMINATED_CLOSING_BRACE,l(),0),v}function Mn(_){return _!==ac&&_!==on}function os(_){C(_),b(_,"'");let v="",T="";for(;v=V(_,Mn);)v==="\\"?T+=L(_):T+=v;const D=_.currentChar();return D===on||D===ca?(i(ds.UNTERMINATED_SINGLE_QUOTE_IN_PLACEHOLDER,l(),0),D===on&&(_.next(),b(_,"'")),T):(b(_,"'"),T)}function L(_){const v=_.currentChar();switch(v){case"\\":case"'":return _.next(),`\\${v}`;case"u":return W(_,v,4);case"U":return W(_,v,6);default:return i(ds.UNKNOWN_ESCAPE_SEQUENCE,l(),0,v),""}}function W(_,v,T){b(_,v);let D="";for(let Z=0;Z<T;Z++){const K=Qs(_);if(!K){i(ds.INVALID_UNICODE_ESCAPE_SEQUENCE,l(),0,`\\${v}${D}${_.currentChar()}`);break}D+=K}return`\\${v}${D}`}function H(_){return _!=="{"&&_!=="}"&&_!==Wn&&_!==on}function J(_){C(_);let v="",T="";for(;v=V(_,H);)T+=v;return T}function hs(_){let v="",T="";for(;v=Q(_);)T+=v;return T}function w(_){const v=T=>{const D=_.currentChar();return D==="{"||D==="@"||D==="|"||D==="("||D===")"||!D||D===Wn?T:(T+=D,_.next(),v(T))};return v("")}function k(_){C(_);const v=b(_,"|");return C(_),v}function A(_,v){let T=null;switch(_.currentChar()){case"{":return v.braceNest>=1&&i(ds.NOT_ALLOW_NEST_PLACEHOLDER,l(),0),_.next(),T=h(v,2,"{"),C(_),v.braceNest++,T;case"}":return v.braceNest>0&&v.currentType===2&&i(ds.EMPTY_PLACEHOLDER,l(),0),_.next(),T=h(v,3,"}"),v.braceNest--,v.braceNest>0&&C(_),v.inLinked&&v.braceNest===0&&(v.inLinked=!1),T;case"@":return v.braceNest>0&&i(ds.UNTERMINATED_CLOSING_BRACE,l(),0),T=F(_,v)||d(v),v.braceNest=0,T;default:{let Z=!0,K=!0,as=!0;if(U(_))return v.braceNest>0&&i(ds.UNTERMINATED_CLOSING_BRACE,l(),0),T=h(v,1,k(_)),v.braceNest=0,v.inLinked=!1,T;if(v.braceNest>0&&(v.currentType===4||v.currentType===5||v.currentType===6))return i(ds.UNTERMINATED_CLOSING_BRACE,l(),0),v.braceNest=0,q(_,v);if(Z=j(_,v))return T=h(v,4,dn(_)),C(_),T;if(K=y(_,v))return T=h(v,5,en(_)),C(_),T;if(as=x(_,v))return T=h(v,6,os(_)),C(_),T;if(!Z&&!K&&!as)return T=h(v,12,J(_)),i(ds.INVALID_TOKEN_IN_PLACEHOLDER,l(),0,T.value),C(_),T;break}}return T}function F(_,v){const{currentType:T}=v;let D=null;const Z=_.currentChar();switch((T===7||T===8||T===11||T===9)&&(Z===on||Z===Wn)&&i(ds.INVALID_LINKED_FORMAT,l(),0),Z){case"@":return _.next(),D=h(v,7,"@"),v.inLinked=!0,D;case".":return C(_),_.next(),h(v,8,".");case":":return C(_),_.next(),h(v,9,":");default:return U(_)?(D=h(v,1,k(_)),v.braceNest=0,v.inLinked=!1,D):R(_,v)||O(_,v)?(C(_),F(_,v)):E(_,v)?(C(_),h(v,11,hs(_))):P(_,v)?(C(_),Z==="{"?A(_,v)||D:h(v,10,w(_))):(T===7&&i(ds.INVALID_LINKED_FORMAT,l(),0),v.braceNest=0,v.inLinked=!1,q(_,v))}}function q(_,v){let T={type:13};if(v.braceNest>0)return A(_,v)||d(v);if(v.inLinked)return F(_,v)||d(v);switch(_.currentChar()){case"{":return A(_,v)||d(v);case"}":return i(ds.UNBALANCED_CLOSING_BRACE,l(),0),_.next(),h(v,3,"}");case"@":return F(_,v)||d(v);default:{if(U(_))return T=h(v,1,k(_)),v.braceNest=0,v.inLinked=!1,T;if(Y(_))return h(v,0,$s(_));break}}return T}function B(){const{currentType:_,offset:v,startLoc:T,endLoc:D}=c;return c.lastType=_,c.lastOffset=v,c.lastStartLoc=T,c.lastEndLoc=D,c.offset=t(),c.startLoc=l(),e.currentChar()===ca?h(c,13):q(e,c)}return{nextToken:B,currentOffset:t,currentPosition:l,context:u}}const tg="parser",lg=/(?:\\\\|\\'|\\u([0-9a-fA-F]{4})|\\U([0-9a-fA-F]{6}))/g;function og(s,n,a){switch(s){case"\\\\":return"\\";case"\\'":return"'";default:{const e=parseInt(n||a,16);return e<=55295||e>=57344?String.fromCodePoint(e):"�"}}}function pg(s={}){const n=s.location!==!1,{onError:a}=s;function e(m,g,j,y,...x){const R=m.currentPosition();if(R.offset+=y,R.column+=y,a){const E=n?Ol(j,R):null,O=ze(g,E,{domain:tg,args:x});a(O)}}function t(m,g,j){const y={type:m};return n&&(y.start=g,y.end=g,y.loc={start:j,end:j}),y}function l(m,g,j,y){n&&(m.end=g,m.loc&&(m.loc.end=j))}function o(m,g){const j=m.context(),y=t(3,j.offset,j.startLoc);return y.value=g,l(y,m.currentOffset(),m.currentPosition()),y}function p(m,g){const j=m.context(),{lastOffset:y,lastStartLoc:x}=j,R=t(5,y,x);return R.index=parseInt(g,10),m.nextToken(),l(R,m.currentOffset(),m.currentPosition()),R}function c(m,g){const j=m.context(),{lastOffset:y,lastStartLoc:x}=j,R=t(4,y,x);return R.key=g,m.nextToken(),l(R,m.currentOffset(),m.currentPosition()),R}function u(m,g){const j=m.context(),{lastOffset:y,lastStartLoc:x}=j,R=t(9,y,x);return R.value=g.replace(lg,og),m.nextToken(),l(R,m.currentOffset(),m.currentPosition()),R}function r(m){const g=m.nextToken(),j=m.context(),{lastOffset:y,lastStartLoc:x}=j,R=t(8,y,x);return g.type!==11?(e(m,ds.UNEXPECTED_EMPTY_LINKED_MODIFIER,j.lastStartLoc,0),R.value="",l(R,y,x),{nextConsumeToken:g,node:R}):(g.value==null&&e(m,ds.UNEXPECTED_LEXICAL_ANALYSIS,j.lastStartLoc,0,In(g)),R.value=g.value||"",l(R,m.currentOffset(),m.currentPosition()),{node:R})}function i(m,g){const j=m.context(),y=t(7,j.offset,j.startLoc);return y.value=g,l(y,m.currentOffset(),m.currentPosition()),y}function h(m){const g=m.context(),j=t(6,g.offset,g.startLoc);let y=m.nextToken();if(y.type===8){const x=r(m);j.modifier=x.node,y=x.nextConsumeToken||m.nextToken()}switch(y.type!==9&&e(m,ds.UNEXPECTED_LEXICAL_ANALYSIS,g.lastStartLoc,0,In(y)),y=m.nextToken(),y.type===2&&(y=m.nextToken()),y.type){case 10:y.value==null&&e(m,ds.UNEXPECTED_LEXICAL_ANALYSIS,g.lastStartLoc,0,In(y)),j.key=i(m,y.value||"");break;case 4:y.value==null&&e(m,ds.UNEXPECTED_LEXICAL_ANALYSIS,g.lastStartLoc,0,In(y)),j.key=c(m,y.value||"");break;case 5:y.value==null&&e(m,ds.UNEXPECTED_LEXICAL_ANALYSIS,g.lastStartLoc,0,In(y)),j.key=p(m,y.value||"");break;case 6:y.value==null&&e(m,ds.UNEXPECTED_LEXICAL_ANALYSIS,g.lastStartLoc,0,In(y)),j.key=u(m,y.value||"");break;default:{e(m,ds.UNEXPECTED_EMPTY_LINKED_KEY,g.lastStartLoc,0);const x=m.context(),R=t(7,x.offset,x.startLoc);return R.value="",l(R,x.offset,x.startLoc),j.key=R,l(j,x.offset,x.startLoc),{nextConsumeToken:y,node:j}}}return l(j,m.currentOffset(),m.currentPosition()),{node:j}}function d(m){const g=m.context(),j=g.currentType===1?m.currentOffset():g.offset,y=g.currentType===1?g.endLoc:g.startLoc,x=t(2,j,y);x.items=[];let R=null;do{const P=R||m.nextToken();switch(R=null,P.type){case 0:P.value==null&&e(m,ds.UNEXPECTED_LEXICAL_ANALYSIS,g.lastStartLoc,0,In(P)),x.items.push(o(m,P.value||""));break;case 5:P.value==null&&e(m,ds.UNEXPECTED_LEXICAL_ANALYSIS,g.lastStartLoc,0,In(P)),x.items.push(p(m,P.value||""));break;case 4:P.value==null&&e(m,ds.UNEXPECTED_LEXICAL_ANALYSIS,g.lastStartLoc,0,In(P)),x.items.push(c(m,P.value||""));break;case 6:P.value==null&&e(m,ds.UNEXPECTED_LEXICAL_ANALYSIS,g.lastStartLoc,0,In(P)),x.items.push(u(m,P.value||""));break;case 7:{const U=h(m);x.items.push(U.node),R=U.nextConsumeToken||null;break}}}while(g.currentType!==13&&g.currentType!==1);const E=g.currentType===1?g.lastOffset:m.currentOffset(),O=g.currentType===1?g.lastEndLoc:m.currentPosition();return l(x,E,O),x}function b(m,g,j,y){const x=m.context();let R=y.items.length===0;const E=t(1,g,j);E.cases=[],E.cases.push(y);do{const O=d(m);R||(R=O.items.length===0),E.cases.push(O)}while(x.currentType!==13);return R&&e(m,ds.MUST_HAVE_MESSAGES_IN_PLURAL,j,0),l(E,m.currentOffset(),m.currentPosition()),E}function f(m){const g=m.context(),{offset:j,startLoc:y}=g,x=d(m);return g.currentType===13?x:b(m,j,y,x)}function C(m){const g=eg(m,Ys({},s)),j=g.context(),y=t(0,j.offset,j.startLoc);return n&&y.loc&&(y.loc.source=m),y.body=f(g),s.onCacheKey&&(y.cacheKey=s.onCacheKey(m)),j.currentType!==13&&e(g,ds.UNEXPECTED_LEXICAL_ANALYSIS,j.lastStartLoc,0,m[j.offset]||""),l(y,g.currentOffset(),g.currentPosition()),y}return{parse:C}}function In(s){if(s.type===13)return"EOF";const n=(s.value||"").replace(/\r?\n/gu,"\\n");return n.length>10?n.slice(0,9)+"…":n}const cg=/<\/?[\w\s="/.':;#-\/]+>/,rg=s=>cg.test(s);function Sn(s){return vs(s)&&xo(s)===0&&(Tn(s,"b")||Tn(s,"body"))}const wi=["b","body"];function ig(s){return wa(s,wi)}const Ci=["c","cases"];function ug(s){return wa(s,Ci,[])}const ki=["s","static"];function hg(s){return wa(s,ki)}const xi=["i","items"];function dg(s){return wa(s,xi,[])}const Si=["t","type"];function xo(s){return wa(s,Si)}const Pi=["v","value"];function Qe(s,n){const a=wa(s,Pi);if(a!=null)return a;throw Ee(n)}const Ei=["m","modifier"];function mg(s){return wa(s,Ei)}const Ti=["k","key"];function jg(s){const n=wa(s,Ti);if(n)return n;throw Ee(6)}function wa(s,n,a){for(let e=0;e<n.length;e++){const t=n[e];if(Tn(s,t)&&s[t]!=null)return s[t]}return a}const Ai=[...wi,...Ci,...ki,...xi,...Ti,...Ei,...Pi,...Si];function Ee(s){return new Error(`unhandled node type: ${s}`)}function il(s){return a=>gg(a,s)}function gg(s,n){const a=ig(n);if(a==null)throw Ee(0);if(xo(a)===1){const l=ug(a);return s.plural(l.reduce((o,p)=>[...o,ec(s,p)],[]))}else return ec(s,a)}function ec(s,n){const a=hg(n);if(a!=null)return s.type==="text"?a:s.normalize([a]);{const e=dg(n).reduce((t,l)=>[...t,Il(s,l)],[]);return s.normalize(e)}}function Il(s,n){const a=xo(n);switch(a){case 3:return Qe(n,a);case 9:return Qe(n,a);case 4:{const e=n;if(Tn(e,"k")&&e.k)return s.interpolate(s.named(e.k));if(Tn(e,"key")&&e.key)return s.interpolate(s.named(e.key));throw Ee(a)}case 5:{const e=n;if(Tn(e,"i")&&Hs(e.i))return s.interpolate(s.list(e.i));if(Tn(e,"index")&&Hs(e.index))return s.interpolate(s.list(e.index));throw Ee(a)}case 6:{const e=n,t=mg(e),l=jg(e);return s.linked(Il(s,l),t?Il(s,t):void 0,s.type)}case 7:return Qe(n,a);case 8:return Qe(n,a);default:throw new Error(`unhandled node on format message part: ${a}`)}}const fg="Detected HTML in '{source}' message. Recommend not using HTML messages to avoid XSS.";function bg(s,n){n&&rg(s)&&va($t(fg,{source:s}))}const _g=s=>s;let Je=Ps();function yg(s,n={}){let a=!1;const e=n.onError||Kj;n.onError=o=>{a=!0,e(o)};const l=pg(n).parse(s);return n.optimize&&Yj(l),n.mangle&&Va(l),{ast:l,detectError:a,code:""}}function vg(s,n){if(ts(s)){const a=Ds(n.warnHtmlMessage)?n.warnHtmlMessage:!0;bg(s,a);const t=(n.onCacheKey||_g)(s),l=Je[t];if(l)return l;const{ast:o,detectError:p}=yg(s,{...n,location:!0,mangle:!1,optimize:!1}),c=il(o);return p?c:Je[t]=c}else{if(!Sn(s))return va(`the message that is resolve with key '${n.key}' is not supported for jit compilation`),(()=>s);const a=s.cacheKey;if(a){const e=Je[a];return e||(Je[a]=il(s))}else return il(s)}}let Te=null;function wg(s){Te=s}function Cg(s,n,a){Te&&Te.emit("i18n:init",{timestamp:Date.now(),i18n:s,version:n,meta:a})}const kg=xg("function:translate");function xg(s){return n=>Te&&Te.emit(s,n)}const sn={INVALID_ARGUMENT:Gj,INVALID_DATE_ARGUMENT:18,INVALID_ISO_DATE_ARGUMENT:19,NOT_SUPPORT_NON_STRING_MESSAGE:20,NOT_SUPPORT_LOCALE_PROMISE_VALUE:21,NOT_SUPPORT_LOCALE_ASYNC_FUNCTION:22,NOT_SUPPORT_LOCALE_TYPE:23},Sg=24;function Zn(s){return ze(s,null,{messages:Pg})}const Pg={[sn.INVALID_ARGUMENT]:"Invalid arguments",[sn.INVALID_DATE_ARGUMENT]:"The date provided is an invalid Date object.Make sure your Date represents a valid date.",[sn.INVALID_ISO_DATE_ARGUMENT]:"The argument provided is not a valid ISO date string",[sn.NOT_SUPPORT_NON_STRING_MESSAGE]:"Not support non-string message",[sn.NOT_SUPPORT_LOCALE_PROMISE_VALUE]:"cannot support promise value",[sn.NOT_SUPPORT_LOCALE_ASYNC_FUNCTION]:"cannot support async function",[sn.NOT_SUPPORT_LOCALE_TYPE]:"cannot support locale type"};function So(s,n){return n.locale!=null?tc(n.locale):tc(s.locale)}let ul;function tc(s){if(ts(s))return s;if(Ls(s)){if(s.resolvedOnce&&ul!=null)return ul;if(s.constructor.name==="Function"){const n=s();if(qj(n))throw Zn(sn.NOT_SUPPORT_LOCALE_PROMISE_VALUE);return ul=n}else throw Zn(sn.NOT_SUPPORT_LOCALE_ASYNC_FUNCTION)}else throw Zn(sn.NOT_SUPPORT_LOCALE_TYPE)}function Eg(s,n,a){return[...new Set([a,...qs(n)?n:vs(n)?Object.keys(n):ts(n)?[n]:[a]])]}function Ri(s,n,a){const e=ts(a)?a:Ct,t=s;t.__localeChainCache||(t.__localeChainCache=new Map);let l=t.__localeChainCache.get(e);if(!l){l=[];let o=[a];for(;qs(o);)o=lc(l,o,n);const p=qs(n)||!ys(n)?n:n.default?n.default:null;o=ts(p)?[p]:p,qs(o)&&lc(l,o,!1),t.__localeChainCache.set(e,l)}return l}function lc(s,n,a){let e=!0;for(let t=0;t<n.length&&Ds(e);t++){const l=n[t];ts(l)&&(e=Tg(s,n[t],a))}return e}function Tg(s,n,a){let e;const t=n.split("-");do{const l=t.join("-");e=Ag(s,l,a),t.splice(-1,1)}while(t.length&&e===!0);return e}function Ag(s,n,a){let e=!1;if(!s.includes(n)&&(e=!0,n)){e=n[n.length-1]!=="!";const t=n.replace(/!/g,"");s.push(t),(qs(a)||ys(a))&&a[t]&&(e=a[t])}return e}const Ca=[];Ca[0]={w:[0],i:[3,0],"[":[4],o:[7]};Ca[1]={w:[1],".":[2],"[":[4],o:[7]};Ca[2]={w:[2],i:[3,0],0:[3,0]};Ca[3]={i:[3,0],0:[3,0],w:[1,1],".":[2,1],"[":[4,1],o:[7,1]};Ca[4]={"'":[5,0],'"':[6,0],"[":[4,2],"]":[1,3],o:8,l:[4,0]};Ca[5]={"'":[4,0],o:8,l:[5,0]};Ca[6]={'"':[4,0],o:8,l:[6,0]};const Rg=/^\s?(?:true|false|-?[\d.]+|'[^']*'|"[^"]*")\s?$/;function Lg(s){return Rg.test(s)}function Mg(s){const n=s.charCodeAt(0),a=s.charCodeAt(s.length-1);return n===a&&(n===34||n===39)?s.slice(1,-1):s}function Dg(s){if(s==null)return"o";switch(s.charCodeAt(0)){case 91:case 93:case 46:case 34:case 39:return s;case 95:case 36:case 45:return"i";case 9:case 10:case 13:case 160:case 65279:case 8232:case 8233:return"w"}return"i"}function Og(s){const n=s.trim();return s.charAt(0)==="0"&&isNaN(parseInt(s))?!1:Lg(n)?Mg(n):"*"+n}function Ig(s){const n=[];let a=-1,e=0,t=0,l,o,p,c,u,r,i;const h=[];h[0]=()=>{o===void 0?o=p:o+=p},h[1]=()=>{o!==void 0&&(n.push(o),o=void 0)},h[2]=()=>{h[0](),t++},h[3]=()=>{if(t>0)t--,e=4,h[0]();else{if(t=0,o===void 0||(o=Og(o),o===!1))return!1;h[1]()}};function d(){const b=s[a+1];if(e===5&&b==="'"||e===6&&b==='"')return a++,p="\\"+b,h[0](),!0}for(;e!==null;)if(a++,l=s[a],!(l==="\\"&&d())){if(c=Dg(l),i=Ca[e],u=i[c]||i.l||8,u===8||(e=u[0],u[1]!==void 0&&(r=h[u[1]],r&&(p=l,r()===!1))))return;if(e===7)return n}}const oc=new Map;function Ng(s,n){return vs(s)?s[n]:null}function Fg(s,n){if(!vs(s))return null;let a=oc.get(n);if(a||(a=Ig(n),a&&oc.set(n,a)),!a)return null;const e=a.length;let t=s,l=0;for(;l<e;){const o=a[l];if(Ai.includes(o)&&Sn(t))return null;const p=t[o];if(p===void 0||Ls(t))return null;t=p,l++}return t}const mn={NOT_FOUND_KEY:1,FALLBACK_TO_TRANSLATE:2,CANNOT_FORMAT_NUMBER:3,FALLBACK_TO_NUMBER_FORMAT:4,CANNOT_FORMAT_DATE:5,FALLBACK_TO_DATE_FORMAT:6,EXPERIMENTAL_CUSTOM_MESSAGE_COMPILER:7},$g=8,Bg={[mn.NOT_FOUND_KEY]:"Not found '{key}' key in '{locale}' locale messages.",[mn.FALLBACK_TO_TRANSLATE]:"Fall back to translate '{key}' key with '{target}' locale.",[mn.CANNOT_FORMAT_NUMBER]:"Cannot format a number value due to not supported Intl.NumberFormat.",[mn.FALLBACK_TO_NUMBER_FORMAT]:"Fall back to number format '{key}' key with '{target}' locale.",[mn.CANNOT_FORMAT_DATE]:"Cannot format a date value due to not supported Intl.DateTimeFormat.",[mn.FALLBACK_TO_DATE_FORMAT]:"Fall back to datetime format '{key}' key with '{target}' locale.",[mn.EXPERIMENTAL_CUSTOM_MESSAGE_COMPILER]:"This project is using Custom Message Compiler, which is an experimental feature. It may receive breaking changes or be removed in the future."};function Oa(s,...n){return $t(Bg[s],...n)}const qg="12.0.0-alpha.3",qt=-1,Ct="en-US",kt="",pc=s=>`${s.charAt(0).toLocaleUpperCase()}${s.substr(1)}`;function zg(){return{upper:(s,n)=>n==="text"&&ts(s)?s.toUpperCase():n==="vnode"&&vs(s)&&"__v_isVNode"in s?s.children.toUpperCase():s,lower:(s,n)=>n==="text"&&ts(s)?s.toLowerCase():n==="vnode"&&vs(s)&&"__v_isVNode"in s?s.children.toLowerCase():s,capitalize:(s,n)=>n==="text"&&ts(s)?pc(s):n==="vnode"&&vs(s)&&"__v_isVNode"in s?pc(s.children):s}}let Li;function Vg(s){Li=s}let Mi;function Ug(s){Mi=s}let Di;function Hg(s){Di=s}let Oi=null;const Gg=s=>{Oi=s},Wg=()=>Oi;let Ii=null;const cc=s=>{Ii=s},Kg=()=>Ii;let rc=0;function Xg(s={}){const n=Ls(s.onWarn)?s.onWarn:va,a=ts(s.version)?s.version:qg,e=ts(s.locale)||Ls(s.locale)?s.locale:Ct,t=Ls(e)?Ct:e,l=qs(s.fallbackLocale)||ys(s.fallbackLocale)||ts(s.fallbackLocale)||s.fallbackLocale===!1?s.fallbackLocale:t,o=ys(s.messages)?s.messages:hl(t),p=ys(s.datetimeFormats)?s.datetimeFormats:hl(t),c=ys(s.numberFormats)?s.numberFormats:hl(t),u=Ys(Ps(),s.modifiers,zg()),r=s.pluralRules||Ps(),i=Ls(s.missing)?s.missing:null,h=Ds(s.missingWarn)||wt(s.missingWarn)?s.missingWarn:!0,d=Ds(s.fallbackWarn)||wt(s.fallbackWarn)?s.fallbackWarn:!0,b=!!s.fallbackFormat,f=!!s.unresolving,C=Ls(s.postTranslation)?s.postTranslation:null,m=ys(s.processor)?s.processor:null,g=Ds(s.warnHtmlMessage)?s.warnHtmlMessage:!0,j=!!s.escapeParameter,y=Ls(s.messageCompiler)?s.messageCompiler:Li;Ls(s.messageCompiler)&&Uj(Oa(mn.EXPERIMENTAL_CUSTOM_MESSAGE_COMPILER));const x=Ls(s.messageResolver)?s.messageResolver:Mi||Ng,R=Ls(s.localeFallbacker)?s.localeFallbacker:Di||Eg,E=vs(s.fallbackContext)?s.fallbackContext:void 0,O=s,P=vs(O.__datetimeFormatters)?O.__datetimeFormatters:new Map,U=vs(O.__numberFormatters)?O.__numberFormatters:new Map,Y=vs(O.__meta)?O.__meta:{};rc++;const V={version:a,cid:rc,locale:e,fallbackLocale:l,messages:o,modifiers:u,pluralRules:r,missing:i,missingWarn:h,fallbackWarn:d,fallbackFormat:b,unresolving:f,postTranslation:C,processor:m,warnHtmlMessage:g,escapeParameter:j,messageCompiler:y,messageResolver:x,localeFallbacker:R,fallbackContext:E,onWarn:n,__meta:Y};return V.datetimeFormats=p,V.numberFormats=c,V.__datetimeFormatters=P,V.__numberFormatters=U,V.__v_emitter=O.__v_emitter!=null?O.__v_emitter:void 0,Cg(V,a,Y),V}const hl=s=>({[s]:Ps()});function zt(s,n){return s instanceof RegExp?s.test(n):s}function Ni(s,n){return s instanceof RegExp?s.test(n):s}function Po(s,n,a,e,t){const{missing:l,onWarn:o}=s;{const p=s.__v_emitter;p&&p.emit("missing",{locale:a,key:n,type:t,groupId:`${t}:${n}`})}if(l!==null){const p=l(s,a,n,t);return ts(p)?p:n}else return Ni(e,n)&&o(Oa(mn.NOT_FOUND_KEY,{key:n,locale:a})),n}function oe(s,n,a){const e=s;e.__localeChainCache=new Map,s.localeFallbacker(s,a,n)}function Fi(s,n){return s===n?!1:s.split("-")[0]===n.split("-")[0]}function Yg(s,n){const a=n.indexOf(s);if(a===-1)return!1;for(let e=a+1;e<n.length;e++)if(Fi(s,n[e]))return!0;return!1}const ic=typeof Intl<"u",$i={dateTimeFormat:ic&&typeof Intl.DateTimeFormat<"u",numberFormat:ic&&typeof Intl.NumberFormat<"u"};function uc(s,...n){const{datetimeFormats:a,unresolving:e,fallbackLocale:t,onWarn:l,localeFallbacker:o}=s,{__datetimeFormatters:p}=s;if(!$i.dateTimeFormat)return l(Oa(mn.CANNOT_FORMAT_DATE)),kt;const[c,u,r,i]=Nl(...n),h=Ds(r.missingWarn)?r.missingWarn:s.missingWarn,d=Ds(r.fallbackWarn)?r.fallbackWarn:s.fallbackWarn,b=!!r.part,f=So(s,r),C=o(s,t,f);if(!ts(c)||c==="")return new Intl.DateTimeFormat(f,i).format(u);let m={},g,j=null,y=f,x=null;const R="datetime format";for(let P=0;P<C.length;P++){if(g=x=C[P],f!==g&&zt(d,c)&&l(Oa(mn.FALLBACK_TO_DATE_FORMAT,{key:c,target:g})),f!==g){const U=s.__v_emitter;U&&U.emit("fallback",{type:R,key:c,from:y,to:x,groupId:`${R}:${c}`})}if(m=a[g]||{},j=m[c],ys(j))break;Po(s,c,g,h,R),y=x}if(!ys(j)||!ts(g))return e?qt:c;let E=`${g}__${c}`;Bt(i)||(E=`${E}__${JSON.stringify(i)}`);let O=p.get(E);return O||(O=new Intl.DateTimeFormat(g,Ys({},j,i)),p.set(E,O)),b?O.formatToParts(u):O.format(u)}const Bi=["localeMatcher","weekday","era","year","month","day","hour","minute","second","timeZoneName","formatMatcher","hour12","timeZone","dateStyle","timeStyle","calendar","dayPeriod","numberingSystem","hourCycle","fractionalSecondDigits"];function Nl(...s){const[n,a,e,t]=s,l=Ps();let o=Ps(),p;if(ts(n)){const c=n.match(/(\d{4}-\d{2}-\d{2})(T|\s)?(.*)/);if(!c)throw Zn(sn.INVALID_ISO_DATE_ARGUMENT);const u=c[3]?c[3].trim().startsWith("T")?`${c[1].trim()}${c[3].trim()}`:`${c[1].trim()}T${c[3].trim()}`:c[1].trim();p=new Date(u);try{p.toISOString()}catch{throw Zn(sn.INVALID_ISO_DATE_ARGUMENT)}}else if(Nj(n)){if(isNaN(n.getTime()))throw Zn(sn.INVALID_DATE_ARGUMENT);p=n}else if(Hs(n))p=n;else throw Zn(sn.INVALID_ARGUMENT);return ts(a)?l.key=a:ys(a)&&Object.keys(a).forEach(c=>{Bi.includes(c)?o[c]=a[c]:l[c]=a[c]}),ts(e)?l.locale=e:ys(e)&&(o=e),ys(t)&&(o=t),[l.key||"",p,l,o]}function hc(s,n,a){const e=s;for(const t in a){const l=`${n}__${t}`;e.__datetimeFormatters.has(l)&&e.__datetimeFormatters.delete(l)}}function dc(s,...n){const{numberFormats:a,unresolving:e,fallbackLocale:t,onWarn:l,localeFallbacker:o}=s,{__numberFormatters:p}=s;if(!$i.numberFormat)return l(Oa(mn.CANNOT_FORMAT_NUMBER)),kt;const[c,u,r,i]=Fl(...n),h=Ds(r.missingWarn)?r.missingWarn:s.missingWarn,d=Ds(r.fallbackWarn)?r.fallbackWarn:s.fallbackWarn,b=!!r.part,f=So(s,r),C=o(s,t,f);if(!ts(c)||c==="")return new Intl.NumberFormat(f,i).format(u);let m={},g,j=null,y=f,x=null;const R="number format";for(let P=0;P<C.length;P++){if(g=x=C[P],f!==g&&zt(d,c)&&l(Oa(mn.FALLBACK_TO_NUMBER_FORMAT,{key:c,target:g})),f!==g){const U=s.__v_emitter;U&&U.emit("fallback",{type:R,key:c,from:y,to:x,groupId:`${R}:${c}`})}if(m=a[g]||{},j=m[c],ys(j))break;Po(s,c,g,h,R),y=x}if(!ys(j)||!ts(g))return e?qt:c;let E=`${g}__${c}`;Bt(i)||(E=`${E}__${JSON.stringify(i)}`);let O=p.get(E);return O||(O=new Intl.NumberFormat(g,Ys({},j,i)),p.set(E,O)),b?O.formatToParts(u):O.format(u)}const qi=["localeMatcher","style","currency","currencyDisplay","currencySign","useGrouping","minimumIntegerDigits","minimumFractionDigits","maximumFractionDigits","minimumSignificantDigits","maximumSignificantDigits","compactDisplay","notation","signDisplay","unit","unitDisplay","roundingMode","roundingPriority","roundingIncrement","trailingZeroDisplay"];function Fl(...s){const[n,a,e,t]=s,l=Ps();let o=Ps();if(!Hs(n))throw Zn(sn.INVALID_ARGUMENT);const p=n;return ts(a)?l.key=a:ys(a)&&Object.keys(a).forEach(c=>{qi.includes(c)?o[c]=a[c]:l[c]=a[c]}),ts(e)?l.locale=e:ys(e)&&(o=e),ys(t)&&(o=t),[l.key||"",p,l,o]}function mc(s,n,a){const e=s;for(const t in a){const l=`${n}__${t}`;e.__numberFormatters.has(l)&&e.__numberFormatters.delete(l)}}const Qg=s=>s,Jg=s=>"",Zg="text",sf=s=>s.length===0?"":yi(s),nf=zj;function jc(s,n){return s=Math.abs(s),n===2?s?s>1?1:0:1:s?Math.min(s,2):0}function af(s){const n=Hs(s.pluralIndex)?s.pluralIndex:-1;return s.named&&(Hs(s.named.count)||Hs(s.named.n))?Hs(s.named.count)?s.named.count:Hs(s.named.n)?s.named.n:n:n}function ef(s,n){n.count||(n.count=s),n.n||(n.n=s)}function tf(s={}){const n=s.locale,a=af(s),e=vs(s.pluralRules)&&ts(n)&&Ls(s.pluralRules[n])?s.pluralRules[n]:jc,t=vs(s.pluralRules)&&ts(n)&&Ls(s.pluralRules[n])?jc:void 0,l=m=>m[e(a,m.length,t)],o=s.list||[],p=m=>o[m],c=s.named||Ps();Hs(s.pluralIndex)&&ef(a,c);const u=m=>c[m];function r(m,g){const j=Ls(s.messages)?s.messages(m,!!g):vs(s.messages)?s.messages[m]:!1;return j||(s.parent?s.parent.message(m):Jg)}const i=m=>s.modifiers?s.modifiers[m]:Qg,h=ys(s.processor)&&Ls(s.processor.normalize)?s.processor.normalize:sf,d=ys(s.processor)&&Ls(s.processor.interpolate)?s.processor.interpolate:nf,b=ys(s.processor)&&ts(s.processor.type)?s.processor.type:Zg,C={list:p,named:u,plural:l,linked:(m,...g)=>{const[j,y]=g;let x="text",R="";g.length===1?vs(j)?(R=j.modifier||R,x=j.type||x):ts(j)&&(R=j||R):g.length===2&&(ts(j)&&(R=j||R),ts(y)&&(x=y||x));const E=r(m,!0)(C),O=x==="vnode"&&qs(E)&&R?E[0]:E;return R?i(R)(O,x):O},message:r,type:b,interpolate:d,normalize:h,values:Ys(Ps(),o,c)};return C}const gc=()=>"",kn=s=>Ls(s);function fc(s,...n){const{fallbackFormat:a,postTranslation:e,unresolving:t,messageCompiler:l,fallbackLocale:o,messages:p}=s,[c,u]=$l(...n),r=Ds(u.missingWarn)?u.missingWarn:s.missingWarn,i=Ds(u.fallbackWarn)?u.fallbackWarn:s.fallbackWarn,h=Ds(u.escapeParameter)?u.escapeParameter:s.escapeParameter,d=!!u.resolvedMessage,b=ts(u.default)||Ds(u.default)?Ds(u.default)?l?c:()=>c:u.default:a?l?c:()=>c:null,f=a||b!=null&&(ts(b)||Ls(b)),C=So(s,u);h&&lf(u);let[m,g,j]=d?[c,C,p[C]||Ps()]:zi(s,c,C,o,i,r),y=m,x=c;if(!d&&!(ts(y)||Sn(y)||kn(y))&&f&&(y=b,x=y),!d&&(!(ts(y)||Sn(y)||kn(y))||!ts(g)))return t?qt:c;if(ts(y)&&s.messageCompiler==null)return va(`The message format compilation is not supported in this build. Because message compiler isn't included. You need to pre-compilation all message format. So translate function return '${c}'.`),c;let R=!1;const E=()=>{R=!0},O=kn(y)?y:Vi(s,c,g,y,x,E);if(R)return y;const P=rf(s,g,j,u),U=tf(P),Y=of(s,O,U),V=e?e(Y,c):Y;{const z={timestamp:Date.now(),key:ts(c)?c:kn(y)?y.key:"",locale:g||(kn(y)?y.locale:""),format:ts(y)?y:kn(y)?y.source:"",message:V};z.meta=Ys({},s.__meta,Wg()||{}),kg(z)}return V}function lf(s){qs(s.list)?s.list=s.list.map(n=>ts(n)?Jp(n):n):vs(s.named)&&Object.keys(s.named).forEach(n=>{ts(s.named[n])&&(s.named[n]=Jp(s.named[n]))})}function zi(s,n,a,e,t,l){const{messages:o,onWarn:p,messageResolver:c,localeFallbacker:u}=s,r=u(s,e,a);let i=Ps(),h,d=null,b=a,f=null;const C="translate";for(let m=0;m<r.length;m++){if(h=f=r[m],a!==h&&!Fi(a,h)&&zt(t,n)&&p(Oa(mn.FALLBACK_TO_TRANSLATE,{key:n,target:h})),a!==h){const x=s.__v_emitter;x&&x.emit("fallback",{type:C,key:n,from:b,to:f,groupId:`${C}:${n}`})}i=o[h]||Ps();let g=null,j,y;if(aa&&(g=window.performance.now(),j="intlify-message-resolve-start",y="intlify-message-resolve-end",wn&&wn(j)),(d=c(i,n))===null&&(d=i[n]),aa){const x=window.performance.now(),R=s.__v_emitter;R&&g&&d&&R.emit("message-resolve",{type:"message-resolve",key:n,message:d,time:x-g,groupId:`${C}:${n}`}),j&&y&&wn&&Da&&(wn(y),Da("intlify message resolve",j,y))}if(ts(d)||Sn(d)||kn(d))break;if(!Yg(h,r)){const x=Po(s,n,h,l,C);x!==n&&(d=x)}b=f}return[d,h,i]}function Vi(s,n,a,e,t,l){const{messageCompiler:o,warnHtmlMessage:p}=s;if(kn(e)){const h=e;return h.locale=h.locale||a,h.key=h.key||n,h}if(o==null){const h=(()=>e);return h.locale=a,h.key=n,h}let c=null,u,r;aa&&(c=window.performance.now(),u="intlify-message-compilation-start",r="intlify-message-compilation-end",wn&&wn(u));const i=o(e,pf(s,a,t,e,p,l));if(aa){const h=window.performance.now(),d=s.__v_emitter;d&&c&&d.emit("message-compilation",{type:"message-compilation",message:e,time:h-c,groupId:`translate:${n}`}),u&&r&&wn&&Da&&(wn(r),Da("intlify message compilation",u,r))}return i.locale=a,i.key=n,i.source=e,i}function of(s,n,a){let e=null,t,l;aa&&(e=window.performance.now(),t="intlify-message-evaluation-start",l="intlify-message-evaluation-end",wn&&wn(t));const o=n(a);if(aa){const p=window.performance.now(),c=s.__v_emitter;c&&e&&c.emit("message-evaluation",{type:"message-evaluation",value:o,time:p-e,groupId:`translate:${n.key}`}),t&&l&&wn&&Da&&(wn(l),Da("intlify message evaluation",t,l))}return o}function $l(...s){const[n,a,e]=s,t=Ps();if(!ts(n)&&!Hs(n)&&!kn(n)&&!Sn(n))throw Zn(sn.INVALID_ARGUMENT);const l=Hs(n)?String(n):(kn(n),n);return Hs(a)?t.plural=a:ts(a)?t.default=a:ys(a)&&!Bt(a)?t.named=a:qs(a)&&(t.list=a),Hs(e)?t.plural=e:ts(e)?t.default=e:ys(e)&&Ys(t,e),[l,t]}function pf(s,n,a,e,t,l){return{locale:n,key:a,warnHtmlMessage:t,onError:o=>{l&&l(o);{const p=cf(e),c=`Message compilation error: ${o.message}`,u=o.location&&p&&Vj(p,o.location.start.offset,o.location.end.offset),r=s.__v_emitter;r&&p&&r.emit("compile-error",{message:p,error:o.message,start:o.location&&o.location.start.offset,end:o.location&&o.location.end.offset,groupId:`translate:${a}`}),console.error(u?`${c}
${u}`:c)}},onCacheKey:o=>Oj(n,a,o)}}function cf(s){if(ts(s))return s;if(s.loc&&s.loc.source)return s.loc.source}function rf(s,n,a,e){const{modifiers:t,pluralRules:l,messageResolver:o,fallbackLocale:p,fallbackWarn:c,missingWarn:u,fallbackContext:r}=s,h={locale:n,modifiers:t,pluralRules:l,messages:(d,b)=>{let f=o(a,d);if(f==null&&(r||b)){const[,,C]=zi(r||s,d,n,p,c,u);f=o(C,d)}if(ts(f)||Sn(f)){let C=!1;const g=Vi(s,d,n,f,d,()=>{C=!0});return C?gc:g}else return kn(f)?f:gc}};return s.processor&&(h.processor=s.processor),e.list&&(h.list=e.list),e.named&&(h.named=e.named),Hs(e.plural)&&(h.pluralIndex=e.plural),h}const uf="12.0.0-alpha.3";function hf(){console.info(`You are running a development build of vue-i18n.
Make sure to use the production build (*.prod.js) when deploying for production.`)}const zs={UNEXPECTED_RETURN_TYPE:Sg,INVALID_ARGUMENT:25,MUST_BE_CALL_SETUP_TOP:26,NOT_INSTALLED:27,REQUIRED_VALUE:28,INVALID_VALUE:29,CANNOT_SETUP_VUE_DEVTOOLS_PLUGIN:30,NOT_INSTALLED_WITH_PROVIDE:31,UNEXPECTED_ERROR:32,DUPLICATE_USE_I18N_CALLING:33};function Bn(s,...n){return ze(s,null,{messages:df,args:n})}const df={[zs.UNEXPECTED_RETURN_TYPE]:"Unexpected return type in composer",[zs.INVALID_ARGUMENT]:"Invalid argument",[zs.MUST_BE_CALL_SETUP_TOP]:"Must be called at the top of a `setup` function",[zs.NOT_INSTALLED]:"Need to install with `app.use` function",[zs.UNEXPECTED_ERROR]:"Unexpected error",[zs.REQUIRED_VALUE]:"Required in value: {0}",[zs.INVALID_VALUE]:"Invalid value",[zs.CANNOT_SETUP_VUE_DEVTOOLS_PLUGIN]:"Cannot setup vue-devtools plugin",[zs.NOT_INSTALLED_WITH_PROVIDE]:"Need to install with `provide` function",[zs.DUPLICATE_USE_I18N_CALLING]:"Duplicate local-scope `useI18n` call detected. Call `useI18n` only once per component."},Bl=Vn("__translateVNode"),ql=Vn("__datetimeParts"),zl=Vn("__numberParts"),Ae=Vn("__enableEmitter"),Vl=Vn("__disableEmitter"),mf=Vn("__setPluralRules"),jf=Vn("__injectWithOption"),Ul=Vn("__dispose"),Wa={FALLBACK_TO_ROOT:$g,NOT_FOUND_PARENT_SCOPE:9,IGNORE_OBJ_FLATTEN:10},gf={[Wa.FALLBACK_TO_ROOT]:"Fall back to {type} '{key}' with root locale.",[Wa.NOT_FOUND_PARENT_SCOPE]:"Not found parent scope. use the global scope.",[Wa.IGNORE_OBJ_FLATTEN]:"Ignore object flatten: '{key}' key has an string value"};function Eo(s,...n){return $t(gf[s],...n)}function Re(s){if(!vs(s)||Sn(s))return s;for(const n in s)if(Tn(s,n))if(!n.includes("."))vs(s[n])&&Re(s[n]);else{const a=n.split("."),e=a.length-1;let t=s,l=!1;for(let o=0;o<e;o++){if(a[o]==="__proto__")throw new Error(`unsafe key: ${a[o]}`);if(a[o]in t||(t[a[o]]=Ps()),!vs(t[a[o]])){va(Eo(Wa.IGNORE_OBJ_FLATTEN,{key:a[o]})),l=!0;break}t=t[a[o]]}if(l||(Sn(t)?Ai.includes(a[e])||delete s[n]:(t[a[e]]=s[n],delete s[n])),!Sn(t)){const o=t[a[e]];vs(o)&&Re(o)}}return s}function Ui(s,n){const{messages:a,__i18n:e,messageResolver:t,flatJson:l}=n,o=ys(a)?a:qs(e)?Ps():{[s]:Ps()};if(qs(e)&&e.forEach(p=>{if("locale"in p&&"resource"in p){const{locale:c,resource:u}=p;c?(o[c]=o[c]||Ps(),pt(u,o[c])):pt(u,o)}else ts(p)&&pt(JSON.parse(p),o)}),t==null&&l)for(const p in o)Tn(o,p)&&Re(o[p]);return o}function Hi(s){return s.type}function ff(s,n,a){let e=vs(n.messages)?n.messages:Ps();"__i18nGlobal"in a&&(e=Ui(s.locale.value,{messages:e,__i18n:a.__i18nGlobal}));const t=Object.keys(e);t.length&&t.forEach(l=>{s.mergeLocaleMessage(l,e[l])});{if(vs(n.datetimeFormats)){const l=Object.keys(n.datetimeFormats);l.length&&l.forEach(o=>{s.mergeDateTimeFormat(o,n.datetimeFormats[o])})}if(vs(n.numberFormats)){const l=Object.keys(n.numberFormats);l.length&&l.forEach(o=>{s.mergeNumberFormat(o,n.numberFormats[o])})}}}function bc(s){return bs(Fe,null,s,0)}const _c="__INTLIFY_META__",yc=()=>[],bf=()=>!1;let vc=0;function wc(s){return((n,a,e,t)=>s(a,e,la()||void 0,t))}const _f=()=>{const s=la();let n=null;return s&&(n=Hi(s)[_c])?{[_c]:n}:null};function Gi(s={}){const{__root:n,__injectWithOption:a}=s,e=n===void 0,t=s.flatJson,l=aa?rs:io;let o=Ds(s.inheritLocale)?s.inheritLocale:!0;const p=l(n&&o?n.locale.value:ts(s.locale)?s.locale:Ct),c=l(n&&o?n.fallbackLocale.value:ts(s.fallbackLocale)||qs(s.fallbackLocale)||ys(s.fallbackLocale)||s.fallbackLocale===!1?s.fallbackLocale:p.value),u=l(Ui(p.value,s)),r=l(ys(s.datetimeFormats)?s.datetimeFormats:{[p.value]:{}}),i=l(ys(s.numberFormats)?s.numberFormats:{[p.value]:{}});let h=n?n.missingWarn:Ds(s.missingWarn)||wt(s.missingWarn)?s.missingWarn:!0,d=n?n.fallbackWarn:Ds(s.fallbackWarn)||wt(s.fallbackWarn)?s.fallbackWarn:!0,b=n?n.fallbackRoot:Ds(s.fallbackRoot)?s.fallbackRoot:!0,f=!!s.fallbackFormat,C=Ls(s.missing)?s.missing:null,m=Ls(s.missing)?wc(s.missing):null,g=Ls(s.postTranslation)?s.postTranslation:null,j=n?n.warnHtmlMessage:Ds(s.warnHtmlMessage)?s.warnHtmlMessage:!0,y=!!s.escapeParameter;const x=n?n.modifiers:ys(s.modifiers)?s.modifiers:{};let R=s.pluralRules||n&&n.pluralRules,E;E=(()=>{e&&cc(null);const M={version:uf,locale:p.value,fallbackLocale:c.value,messages:u.value,modifiers:x,pluralRules:R,missing:m===null?void 0:m,missingWarn:h,fallbackWarn:d,fallbackFormat:f,unresolving:!0,postTranslation:g===null?void 0:g,warnHtmlMessage:j,escapeParameter:y,messageResolver:s.messageResolver,messageCompiler:s.messageCompiler,__meta:{framework:"vue"}};M.datetimeFormats=r.value,M.numberFormats=i.value,M.__datetimeFormatters=ys(E)?E.__datetimeFormatters:void 0,M.__numberFormatters=ys(E)?E.__numberFormatters:void 0,M.__v_emitter=ys(E)?E.__v_emitter:void 0;const $=Xg(M);return e&&cc($),$})(),oe(E,p.value,c.value);function P(){return[p.value,c.value,u.value,r.value,i.value]}const U=ns({get:()=>p.value,set:M=>{E.locale=M,p.value=M}}),Y=ns({get:()=>c.value,set:M=>{E.fallbackLocale=M,c.value=M,oe(E,p.value,M)}}),V=ns(()=>u.value),z=ns(()=>Object.keys(u.value).sort()),Q=ns(()=>r.value),cs=ns(()=>i.value);function ls(){return Ls(g)?g:null}function X(M){g=M,E.postTranslation=M}function us(){return C}function Fs(M){M!==null&&(m=wc(M)),C=M,E.missing=m}function Qs(M,$){return M!=="translate"||!$.resolvedMessage}const Ms=(M,$,ps,_s,Bs,tn)=>{P();let Ws;try{e||(E.fallbackContext=n?Kg():void 0),Ws=M(E)}finally{e||(E.fallbackContext=void 0)}if(ps!=="translate exists"&&Hs(Ws)&&Ws===qt||ps==="translate exists"&&!Ws){const[fn,Ve]=$();if(n&&ts(fn)&&Qs(ps,Ve)){b&&(zt(d,fn)||Ni(h,fn))&&va(Eo(Wa.FALLBACK_TO_ROOT,{key:fn,type:ps}));{const{__v_emitter:Js}=E;Js&&b&&Js.emit("fallback",{type:ps,key:fn,to:"global",groupId:`${ps}:${fn}`})}}return n&&b?_s(n):Bs(fn)}else{if(tn(Ws))return Ws;throw Bn(zs.UNEXPECTED_RETURN_TYPE)}};function $s(...M){return Ms($=>Reflect.apply(fc,null,[$,...M]),()=>$l(...M),"translate",$=>Reflect.apply($.t,$,[...M]),$=>$,$=>ts($))}function dn(...M){const[$,ps,_s]=M;if(_s&&!vs(_s))throw Bn(zs.INVALID_ARGUMENT);return $s($,ps,Ys({resolvedMessage:!0},_s||{}))}function en(...M){return Ms($=>Reflect.apply(uc,null,[$,...M]),()=>Nl(...M),"datetime format",$=>Reflect.apply($.d,$,[...M]),()=>kt,$=>ts($)||qs($))}function Mn(...M){return Ms($=>Reflect.apply(dc,null,[$,...M]),()=>Fl(...M),"number format",$=>Reflect.apply($.n,$,[...M]),()=>kt,$=>ts($)||qs($))}function os(M){return M.map($=>ts($)||Hs($)||Ds($)?bc(String($)):$)}const W={normalize:os,interpolate:M=>M,type:"vnode"};function H(...M){return Ms($=>{let ps;const _s=$;try{_s.processor=W,ps=Reflect.apply(fc,null,[_s,...M])}finally{_s.processor=null}return ps},()=>$l(...M),"translate",$=>$[Bl](...M),$=>[bc($)],$=>qs($))}function J(...M){return Ms($=>Reflect.apply(dc,null,[$,...M]),()=>Fl(...M),"number format",$=>$[zl](...M),yc,$=>ts($)||qs($))}function hs(...M){return Ms($=>Reflect.apply(uc,null,[$,...M]),()=>Nl(...M),"datetime format",$=>$[ql](...M),yc,$=>ts($)||qs($))}function w(M){R=M,E.pluralRules=R}function k(M,$){return Ms(()=>{if(!M)return!1;const ps=ts($)?$:p.value,_s=q(ps),Bs=E.messageResolver(_s,M);return Sn(Bs)||kn(Bs)||ts(Bs)},()=>[M],"translate exists",ps=>Reflect.apply(ps.te,ps,[M,$]),bf,ps=>Ds(ps))}function A(M){let $=null;const ps=Ri(E,c.value,p.value);for(let _s=0;_s<ps.length;_s++){const Bs=u.value[ps[_s]]||{},tn=E.messageResolver(Bs,M);if(tn!=null){$=tn;break}}return $}function F(M){const $=A(M);return $??(n?n.tm(M)||{}:{})}function q(M){return u.value[M]||{}}function B(M,$){if(t){const ps={[M]:$};for(const _s in ps)Tn(ps,_s)&&Re(ps[_s]);$=ps[M]}u.value[M]=$,E.messages=u.value}function _(M,$){u.value[M]=u.value[M]||{};const ps={[M]:$};if(t)for(const _s in ps)Tn(ps,_s)&&Re(ps[_s]);$=ps[M],pt($,u.value[M]),E.messages=u.value}function v(M){return r.value[M]||{}}function T(M,$){r.value[M]=$,E.datetimeFormats=r.value,hc(E,M,$)}function D(M,$){r.value[M]=Ys(r.value[M]||{},$),E.datetimeFormats=r.value,hc(E,M,$)}function Z(M){return i.value[M]||{}}function K(M,$){i.value[M]=$,E.numberFormats=i.value,mc(E,M,$)}function as(M,$){i.value[M]=Ys(i.value[M]||{},$),E.numberFormats=i.value,mc(E,M,$)}vc++,n&&aa&&(Vs(n.locale,M=>{o&&(p.value=M,E.locale=M,oe(E,p.value,c.value))}),Vs(n.fallbackLocale,M=>{o&&(c.value=M,E.fallbackLocale=M,oe(E,p.value,c.value))}));const es={id:vc,locale:U,fallbackLocale:Y,get inheritLocale(){return o},set inheritLocale(M){o=M,M&&n&&(p.value=n.locale.value,c.value=n.fallbackLocale.value,oe(E,p.value,c.value))},availableLocales:z,messages:V,get modifiers(){return x},get pluralRules(){return R||{}},get isGlobal(){return e},get missingWarn(){return h},set missingWarn(M){h=M,E.missingWarn=h},get fallbackWarn(){return d},set fallbackWarn(M){d=M,E.fallbackWarn=d},get fallbackRoot(){return b},set fallbackRoot(M){b=M},get fallbackFormat(){return f},set fallbackFormat(M){f=M,E.fallbackFormat=f},get warnHtmlMessage(){return j},set warnHtmlMessage(M){j=M,E.warnHtmlMessage=M},get escapeParameter(){return y},set escapeParameter(M){y=M,E.escapeParameter=M},t:$s,getLocaleMessage:q,setLocaleMessage:B,mergeLocaleMessage:_,getPostTranslationHandler:ls,setPostTranslationHandler:X,getMissingHandler:us,setMissingHandler:Fs,[mf]:w};return es.datetimeFormats=Q,es.numberFormats=cs,es.rt=dn,es.te=k,es.tm=F,es.d=en,es.n=Mn,es.getDateTimeFormat=v,es.setDateTimeFormat=T,es.mergeDateTimeFormat=D,es.getNumberFormat=Z,es.setNumberFormat=K,es.mergeNumberFormat=as,es[jf]=a,es[Bl]=H,es[ql]=hs,es[zl]=J,es[Ae]=M=>{E.__v_emitter=M},es[Vl]=()=>{E.__v_emitter=void 0},es}var _e=typeof global<"u"?global:typeof self<"u"?self:typeof window<"u"?window:{};function yf(){return Wi().__VUE_DEVTOOLS_GLOBAL_HOOK__}function Wi(){return typeof navigator<"u"&&typeof window<"u"?window:typeof _e<"u"?_e:{}}const vf=typeof Proxy=="function",wf="devtools-plugin:setup",Cf="plugin:settings:set";let Fa,Hl;function kf(){var s;return Fa!==void 0||(typeof window<"u"&&window.performance?(Fa=!0,Hl=window.performance):typeof _e<"u"&&(!((s=_e.perf_hooks)===null||s===void 0)&&s.performance)?(Fa=!0,Hl=_e.perf_hooks.performance):Fa=!1),Fa}function xf(){return kf()?Hl.now():Date.now()}class Sf{constructor(n,a){this.target=null,this.targetQueue=[],this.onQueue=[],this.plugin=n,this.hook=a;const e={};if(n.settings)for(const o in n.settings){const p=n.settings[o];e[o]=p.defaultValue}const t=`__vue-devtools-plugin-settings__${n.id}`;let l=Object.assign({},e);try{const o=localStorage.getItem(t),p=JSON.parse(o);Object.assign(l,p)}catch{}this.fallbacks={getSettings(){return l},setSettings(o){try{localStorage.setItem(t,JSON.stringify(o))}catch{}l=o},now(){return xf()}},a&&a.on(Cf,(o,p)=>{o===this.plugin.id&&this.fallbacks.setSettings(p)}),this.proxiedOn=new Proxy({},{get:(o,p)=>this.target?this.target.on[p]:(...c)=>{this.onQueue.push({method:p,args:c})}}),this.proxiedTarget=new Proxy({},{get:(o,p)=>this.target?this.target[p]:p==="on"?this.proxiedOn:Object.keys(this.fallbacks).includes(p)?(...c)=>(this.targetQueue.push({method:p,args:c,resolve:()=>{}}),this.fallbacks[p](...c)):(...c)=>new Promise(u=>{this.targetQueue.push({method:p,args:c,resolve:u})})})}async setRealTarget(n){this.target=n;for(const a of this.onQueue)this.target.on[a.method](...a.args);for(const a of this.targetQueue)a.resolve(await this.target[a.method](...a.args))}}function Pf(s,n){const a=s,e=Wi(),t=yf(),l=vf&&a.enableEarlyProxy;if(t&&(e.__VUE_DEVTOOLS_PLUGIN_API_AVAILABLE__||!l))t.emit(wf,s,n);else{const o=l?new Sf(a,t):null;(e.__VUE_DEVTOOLS_PLUGINS__=e.__VUE_DEVTOOLS_PLUGINS__||[]).push({pluginDescriptor:a,setupFn:n,proxy:o}),o&&n(o.proxiedTarget)}}const Ki="vue-i18n: composer properties",dl={"vue-devtools-plugin-vue-i18n":"Vue I18n DevTools","vue-i18n-resource-inspector":"Vue I18n DevTools","vue-i18n-timeline":"Vue I18n"},Ef={"vue-i18n-resource-inspector":"Search for scopes ..."},Tf={"vue-i18n-timeline":16764185};let Gl;async function Af(s,n){return new Promise((a,e)=>{try{Pf({id:"vue-devtools-plugin-vue-i18n",label:dl["vue-devtools-plugin-vue-i18n"],packageName:"vue-i18n",homepage:"https://vue-i18n.intlify.dev",logo:"https://vue-i18n.intlify.dev/vue-i18n-devtools-logo.png",componentStateTypes:[Ki],app:s},t=>{Gl=t,t.on.visitComponentTree(({componentInstance:o,treeNode:p})=>{Rf(o,p,n)}),t.on.inspectComponent(({componentInstance:o,instanceData:p})=>{o.vnode.el&&o.vnode.el.__VUE_I18N__&&p&&Lf(p,o.vnode.el.__VUE_I18N__)}),t.addInspector({id:"vue-i18n-resource-inspector",label:dl["vue-i18n-resource-inspector"],icon:"language",treeFilterPlaceholder:Ef["vue-i18n-resource-inspector"]}),t.on.getInspectorTree(o=>{o.app===s&&o.inspectorId==="vue-i18n-resource-inspector"&&Nf(o,n)});const l=new Map;t.on.getInspectorState(async o=>{if(o.app===s&&o.inspectorId==="vue-i18n-resource-inspector")if(t.unhighlightElement(),$f(o,n),o.nodeId==="global"){if(!l.has(o.app)){const[p]=await t.getComponentInstances(o.app);l.set(o.app,p)}t.highlightElement(l.get(o.app))}else{const p=Ff(o.nodeId,n);p&&t.highlightElement(p)}}),t.on.editInspectorState(o=>{o.app===s&&o.inspectorId==="vue-i18n-resource-inspector"&&qf(o,n)}),t.addTimelineLayer({id:"vue-i18n-timeline",label:dl["vue-i18n-timeline"],color:Tf["vue-i18n-timeline"]}),a(!0)})}catch(t){console.error(t),e(!1)}})}function Xi(s){return s.type.name||s.type.displayName||s.type.__file||"Anonymous"}function Rf(s,n,a){const e=a.global;if(s&&s.vnode.el&&s.vnode.el.__VUE_I18N__&&s.vnode.el.__VUE_I18N__!==e){const t={label:`i18n (${Xi(s)} Scope)`,textColor:0,backgroundColor:16764185};n.tags.push(t)}}function Lf(s,n){const a=Ki;s.state.push({type:a,key:"locale",editable:!0,value:n.locale.value}),s.state.push({type:a,key:"availableLocales",editable:!1,value:n.availableLocales}),s.state.push({type:a,key:"fallbackLocale",editable:!0,value:n.fallbackLocale.value}),s.state.push({type:a,key:"inheritLocale",editable:!0,value:n.inheritLocale}),s.state.push({type:a,key:"messages",editable:!1,value:To(n.messages.value)}),s.state.push({type:a,key:"datetimeFormats",editable:!1,value:n.datetimeFormats.value}),s.state.push({type:a,key:"numberFormats",editable:!1,value:n.numberFormats.value})}function To(s){const n={};return Object.keys(s).forEach(a=>{const e=s[a];Ls(e)&&"source"in e?n[a]=If(e):Sn(e)&&e.loc&&e.loc.source?n[a]=e.loc.source:vs(e)?n[a]=To(e):n[a]=e}),n}const Mf={"<":"&lt;",">":"&gt;",'"':"&quot;","&":"&amp;"};function Df(s){return s.replace(/[<>"&]/g,Of)}function Of(s){return Mf[s]||s}function If(s){return{_custom:{type:"function",display:`<span>ƒ</span> ${s.source?`("${Df(s.source)}")`:"(?)"}`}}}function Nf(s,n){s.rootNodes.push({id:"global",label:"Global Scope"});const a=n.global;for(const[e,t]of n.__instances){const l=t;a!==l&&s.rootNodes.push({id:l.id.toString(),label:`${Xi(e)} Scope`})}}function Ff(s,n){let a=null;if(s!=="global"){for(const[e,t]of n.__instances.entries())if(t.id.toString()===s){a=e;break}}return a}function Yi(s,n){if(s==="global")return n.global;{const a=Array.from(n.__instances.values()).find(e=>e.id.toString()===s);return a||null}}function $f(s,n){const a=Yi(s.nodeId,n);return a&&(s.state=Bf(a)),null}function Bf(s){const n={},a="Locale related info",e=[{type:a,key:"locale",editable:!0,value:s.locale.value},{type:a,key:"fallbackLocale",editable:!0,value:s.fallbackLocale.value},{type:a,key:"availableLocales",editable:!1,value:s.availableLocales},{type:a,key:"inheritLocale",editable:!0,value:s.inheritLocale}];n[a]=e;const t="Locale messages info",l=[{type:t,key:"messages",editable:!1,value:To(s.messages.value)}];n[t]=l;{const o="Datetime formats info",p=[{type:o,key:"datetimeFormats",editable:!1,value:s.datetimeFormats.value}];n[o]=p;const c="Datetime formats info",u=[{type:c,key:"numberFormats",editable:!1,value:s.numberFormats.value}];n[c]=u}return n}function Wl(s,n){if(Gl){let a;n&&"groupId"in n&&(a=n.groupId,delete n.groupId),Gl.addTimelineEvent({layerId:"vue-i18n-timeline",event:{title:s,groupId:a,time:Date.now(),meta:{},data:n||{},logType:s==="compile-error"?"error":s==="fallback"||s==="missing"?"warning":"default"}})}}function qf(s,n){const a=Yi(s.nodeId,n);if(a){const[e]=s.path;e==="locale"&&ts(s.state.value)?a.locale.value=s.state.value:e==="fallbackLocale"&&(ts(s.state.value)||qs(s.state.value)||vs(s.state.value))?a.fallbackLocale.value=s.state.value:e==="inheritLocale"&&Ds(s.state.value)&&(a.inheritLocale=s.state.value)}}const Ao={tag:{type:[String,Object]},locale:{type:String},scope:{type:String,validator:s=>s==="parent"||s==="global",default:"parent"},i18n:{type:Object}};function zf({slots:s},n){return n.length===1&&n[0]==="default"?(s.default?s.default():[]).reduce((e,t)=>[...e,...t.type===gs?t.children:[t]],[]):n.reduce((a,e)=>{const t=s[e];return t&&(a[e]=t()),a},Ps())}function Qi(){return gs}function Vf(s){return qs(s)&&!ts(s[0])}function Ji(s,n,a,e){const{slots:t,attrs:l}=n;return()=>{const o={part:!0};let p=Ps();s.locale&&(o.locale=s.locale),ts(s.format)?o.key=s.format:vs(s.format)&&(ts(s.format.key)&&(o.key=s.format.key),p=Object.keys(s.format).reduce((h,d)=>a.includes(d)?Ys(Ps(),h,{[d]:s.format[d]}):h,Ps()));const c=e(s.value,o,p);let u=[o.key];qs(c)?u=c.map((h,d)=>{const b=t[h.type],f=b?b({[h.type]:h.value,index:d,parts:c}):[h.value];return Vf(f)&&(f[0].key=`${h.type}-${d}`),f}):ts(c)&&(u=[c]);const r=Ys(Ps(),l),i=ts(s.tag)||vs(s.tag)?s.tag:Qi();return Be(i,r,u)}}const Uf=Ns({name:"i18n-d",props:Ys({value:{type:[Number,Date],required:!0},format:{type:[String,Object]}},Ao),setup(s,n){const a=s.i18n||Us({useScope:s.scope,__useComponent:!0});return Ji(s,n,Bi,(...e)=>a[ql](...e))}}),Cc=Uf,Hf=Ns({name:"i18n-n",props:Ys({value:{type:Number,required:!0},format:{type:[String,Object]}},Ao),setup(s,n){const a=s.i18n||Us({useScope:s.scope,__useComponent:!0});return Ji(s,n,qi,(...e)=>a[zl](...e))}}),kc=Hf,Gf=Ns({name:"i18n-t",props:Ys({},{keypath:{type:String,required:!0},plural:{type:[Number,String],validator:s=>Hs(s)||!isNaN(s)}},Ao),setup(s,n){const{slots:a,attrs:e}=n,t=s.i18n||Us({useScope:s.scope,__useComponent:!0});return()=>{const l=Object.keys(a).filter(i=>i[0]!=="_"),o=Ps();s.locale&&(o.locale=s.locale),s.plural!==void 0&&(o.plural=ts(s.plural)?+s.plural:s.plural);const p=zf(n,l),c=t[Bl](s.keypath,p,o),u=Ys(Ps(),e),r=ts(s.tag)||vs(s.tag)?s.tag:Qi();return Be(r,u,c)}}}),xc=Gf;function Wf(s,...n){const a=ys(n[0])?n[0]:{};(Ds(a.globalInstall)?a.globalInstall:!0)&&([xc.name,"I18nT"].forEach(t=>s.component(t,xc)),[kc.name,"I18nN"].forEach(t=>s.component(t,kc)),[Cc.name,"I18nD"].forEach(t=>s.component(t,Cc)))}const Kf=Vn("global-vue-i18n");function Xf(s={}){const n=Ds(s.globalInjection)?s.globalInjection:!0,a=new Map,[e,t]=Yf(s),l=Vn("vue-i18n");function o(r){return a.get(r)||null}function p(r,i){a.set(r,i)}function c(r){a.delete(r)}const u={async install(r,...i){if(r.__VUE_I18N__=u,r.__VUE_I18N_SYMBOL__=l,r.provide(r.__VUE_I18N_SYMBOL__,u),ys(i[0])){const b=i[0];u.__composerExtend=b.__composerExtend}let h=null;n&&(h=tb(r,u.global)),Wf(r,...i);const d=r.unmount;r.unmount=()=>{h&&h(),u.dispose(),d()};{if(!await Af(r,u))throw Bn(zs.CANNOT_SETUP_VUE_DEVTOOLS_PLUGIN);const f=vi(),C=t;C[Ae]&&C[Ae](f),f.on("*",Wl)}},get global(){return t},dispose(){e.stop()},__instances:a,__getInstance:o,__setInstance:p,__deleteInstance:c};return u}function Us(s={}){const n=la();if(n==null)throw Bn(zs.MUST_BE_CALL_SETUP_TOP);if(!n.isCE&&n.appContext.app!=null&&!n.appContext.app.__VUE_I18N_SYMBOL__)throw Bn(zs.NOT_INSTALLED);const a=Qf(n),e=Zf(a),t=Hi(n),l=Jf(s,t);if(l==="global")return ff(e,s,t),e;if(l==="parent"){let c=sb(a,n,s.__useComponent);return c==null&&(va(Eo(Wa.NOT_FOUND_PARENT_SCOPE)),c=e),c}const o=a;let p=o.__getInstance(n);if(p==null){const c=Ys({},s);"__i18n"in t&&(c.__i18n=t.__i18n),e&&(c.__root=e),p=Gi(c),o.__composerExtend&&(p[Ul]=o.__composerExtend(p)),ab(o,n,p),o.__setInstance(n,p)}else if(l==="local")throw Bn(zs.DUPLICATE_USE_I18N_CALLING);return p}function Yf(s){const n=no(),a=n.run(()=>Gi(s));if(a==null)throw Bn(zs.UNEXPECTED_ERROR);return[n,a]}function Qf(s){const n=gn(s.isCE?Kf:s.appContext.app.__VUE_I18N_SYMBOL__);if(!n)throw Bn(s.isCE?zs.NOT_INSTALLED_WITH_PROVIDE:zs.UNEXPECTED_ERROR);return n}function Jf(s,n){return Bt(s)?"__i18n"in n?"local":"global":s.useScope?s.useScope:"local"}function Zf(s){return s.global}function sb(s,n,a=!1){let e=null;const t=n.root;let l=nb(n,a);for(;l!=null&&(e=s.__getInstance(l),!(e!=null||t===l));)l=l.parent;return e}function nb(s,n=!1){return s==null?null:n&&s.vnode.ctx||s.parent}function ab(s,n,a){let e=null;En(()=>{if(n.vnode.el){n.vnode.el.__VUE_I18N__=a,e=vi();const t=a;t[Ae]&&t[Ae](e),e.on("*",Wl)}},n),mo(()=>{const t=a;n.vnode.el&&n.vnode.el.__VUE_I18N__&&(e&&e.off("*",Wl),t[Vl]&&t[Vl](),delete n.vnode.el.__VUE_I18N__),s.__deleteInstance(n);const l=t[Ul];l&&(l(),delete t[Ul])},n)}const eb=["locale","fallbackLocale","availableLocales"],Sc=["t","rt","d","n","tm","te"];function tb(s,n){const a=Object.create(null);return eb.forEach(t=>{const l=Object.getOwnPropertyDescriptor(n,t);if(!l)throw Bn(zs.UNEXPECTED_ERROR);const o=Os(l.value)?{get(){return l.value.value},set(p){l.value.value=p}}:{get(){return l.get&&l.get()}};Object.defineProperty(a,t,o)}),s.config.globalProperties.$i18n=a,Sc.forEach(t=>{const l=Object.getOwnPropertyDescriptor(n,t);if(!l||!l.value)throw Bn(zs.UNEXPECTED_ERROR);Object.defineProperty(s.config.globalProperties,`$${t}`,l)}),()=>{delete s.config.globalProperties.$i18n,Sc.forEach(t=>{delete s.config.globalProperties[`$${t}`]})}}Vg(vg);Ug(Fg);Hg(Ri);{const s=$j();s.__INTLIFY__=!0,wg(s.__INTLIFY_DEVTOOLS_GLOBAL_HOOK__)}hf();const lb={tags:"标签",articles:"文章",words:"字数",prevPage:"上一页",nextPage:"下一页",pagination:"分页导航",search:"搜索",theme:"主题",language:"语言",menu:"菜单",close:"关闭",searchPlaceholder:"搜索文章标题或正文...",searchNoResults:"未找到匹配的文章",searchUnavailable:"搜索索引加载失败",searchEscHint:"Esc 关闭",searchResultsLabel:"条结果",categories:"分类",resources:"资源",about:"关于",greeting:"你好，我是",greetingPrefix:"//",developer:"开发者",wordUnit:"字",recentPosts:"最近文章",notes:"笔记",projects:"项目",topics:"课题",seeMore:"查看更多",countPosts:"{count} 篇",countProjects:"{count} 个项目",countTopics:"{count} 个课题",countWords:"{count} 字",copyCode:"复制代码",copyFailed:"复制失败",copyArticle:"复制文章",copyTable:"复制表格",copied:"已复制",anchorHeading:"置顶当前标题",uncategorized:"未分类",updatedAt:"更新于",postReadingTime:"{minutes} 分钟",articleReadingTime:"阅读约 {minutes} 分钟",tableOfContents:"本页目录",openToc:"打开目录",backToTop:"回到顶部",resourceSubtitle:"生物信息学与结构生物学领域常用工具",pageNotFound:"页面不存在或已被移动",backHome:"返回首页",experience:"经历",thanks:"感谢您的关注！",designedByPrefix:"由",designedBySuffix:"设计",clearFilter:"清除筛选",backToArticle:"返回文章"},ob={tags:"Tags",articles:"Articles",words:"Words",prevPage:"Previous",nextPage:"Next",pagination:"Pagination",search:"Search",theme:"Theme",language:"Language",menu:"Menu",close:"Close",searchPlaceholder:"Search articles by title or content...",searchNoResults:"No matching articles found",searchUnavailable:"Search index failed to load",searchEscHint:"Esc to close",searchResultsLabel:"results",categories:"Categories",resources:"Resources",about:"About",greeting:"Hello, I'm",greetingPrefix:"//",developer:"Developer",wordUnit:"words",recentPosts:"Recent Posts",notes:"Notes",projects:"Projects",topics:"Topics",seeMore:"See More",countPosts:"{count} posts",countProjects:"{count} projects",countTopics:"{count} topics",countWords:"{count} words",copyCode:"Copy code",copyFailed:"Copy failed",copyArticle:"Copy article",copyTable:"Copy table",copied:"Copied",anchorHeading:"Anchor to heading",uncategorized:"Uncategorized",updatedAt:"Updated at",postReadingTime:"{minutes} min",articleReadingTime:"Reading about {minutes} minutes",tableOfContents:"On this page",openToc:"Open table of contents",backToTop:"Back to top",resourceSubtitle:"Commonly used tools in bioinformatics and structural biology",pageNotFound:"This page could not be found",backHome:"Back to home",experience:"Experience",thanks:"Thank you for your attention!",designedByPrefix:"Designed by",designedBySuffix:"",clearFilter:"Clear filter",backToArticle:"Back to article"},cn={url:"https://zorrooz.github.io",author:"zorrooz",email:"zorrooz@163.com",github:"https://github.com/zorrooz",startYear:2025},Ra=["zh-CN","en-US"],pb=["zh","en"],ml=["auto","light","dark"],Ro="zh-CN",cb={zh:"zh-CN",en:"en-US"},Kl={"zh-CN":"zh","en-US":"en"},xt={"zh-CN":"zh-CN","en-US":"en"},Vt=s=>Kl[s],Lo=s=>{const n=s.match(/^\/(zh|en)(?=\/|$)/);return n?cb[n[1]]:null},Zi=s=>s.replace(/^\/(zh|en)(?=\/|$)/,""),ea=s=>s===Ra[1]?Ra[1]:Ra[0],$a=()=>{if(typeof window<"u"){const s=localStorage.getItem("locale");if(s===Ra[1])return Kl[Ra[1]];if(s===Ra[0])return Kl[Ra[0]];if(navigator.language&&navigator.language.toLowerCase().startsWith("en"))return"en"}return"zh"},su=typeof window<"u"?ea(localStorage.getItem("locale")):Ro;typeof document<"u"&&(document.documentElement.lang=xt[su]);const Le=Xf({locale:su,fallbackLocale:Ro,messages:{"zh-CN":lb,"en-US":ob},globalInjection:!0}),Mo=Zr("locale",()=>{const s=rs(typeof window<"u"?ea(localStorage.getItem("locale")):Ro);return{locale:s,setLocale:e=>{s.value=e,Le.global.locale.value=e,typeof window<"u"&&localStorage.setItem("locale",e),typeof document<"u"&&(document.documentElement.lang=xt[e])},initLocale:()=>{if(typeof window<"u"){const e=Lo(window.location.pathname);e&&(s.value=e)}Le.global.locale.value=s.value,typeof document<"u"&&(document.documentElement.lang=xt[s.value])}}}),Pc=s=>{if(typeof window>"u"||typeof document>"u")return;const n=document.documentElement;if(s==="auto"){const a=window.matchMedia("(prefers-color-scheme: dark)").matches;n.setAttribute("data-bs-theme",a?"dark":"light"),localStorage.removeItem("theme")}else n.setAttribute("data-bs-theme",s),localStorage.setItem("theme",s)},Do=Zr("theme",()=>{const s=rs(typeof window<"u"&&localStorage.getItem("theme")||"auto");return{theme:s,toggleTheme:()=>{const t=(ml.indexOf(s.value)+1)%ml.length;s.value=ml[t],Pc(s.value)},initTheme:()=>{Pc(s.value)}}}),Qa=(s="smooth")=>{typeof window>"u"||window.scrollTo({top:0,behavior:s})};function Ja(s,n){const a=document.documentElement;a&&(a.style.overflow=s,a.style.overscrollBehavior=n)}function Ut(s,n){const a=document.body;a&&(a.style.overflow=s,a.style.overscrollBehavior=n)}function rb(){Ja("hidden","contain"),Ut("hidden","contain")}function Ec(){Ja("",""),Ut("","")}function ib(){const s=window.scrollY||window.pageYOffset||document.documentElement.scrollTop||0;try{const n=document.body;n&&(n.style.position="fixed",n.style.top=`-${s}px`,n.style.left="0",n.style.right="0",n.style.overflow="hidden"),Ja("","contain")}catch{Ja("hidden","contain"),Ut("hidden","contain")}return s}function Tc(s){try{const n=document.body;n&&(n.style.position="",n.style.top="",n.style.left="",n.style.right="",n.style.overflow=""),Ja("",""),typeof s=="number"&&window.scrollTo(0,s)}catch{Ja("",""),Ut("","")}}const Xl=s=>{const n=s.replace(/^\/+/,"").split("/"),a=n[0]==="article"?1:0;return`${n[a]}/${n[a+1]}/${n.slice(a+2).join("/")}.md`},nu=s=>`/article/${s.replace(/\.md$/,"")}`,ct=s=>Array.isArray(s)?s.join("/"):typeof s=="string"&&s.length?s:"",Jn=s=>`/${Vt(ea(Le.global.locale.value))}${s.startsWith("/")?s:`/${s}`}`,ub=(s,n)=>{const e=(Lo(n.path)??ea(Le.global.locale.value))==="zh-CN"?"en-US":"zh-CN";s.push({path:`/${Vt(e)}${Zi(n.path)}`,query:n.query})};function Oo(s,n,a){n&&(s.push({path:`/${Vt(a.locale)}/`,query:{...a.query??{},tag:n}}).catch(()=>{}),a.scroll!==!1&&bn(()=>Qa()))}function fa(s,n){const a=rs(n),{locale:e}=Us(),t=()=>{try{a.value=s()}catch(l){console.error("Failed to load localized content:",l),a.value=n}};return t(),Vs(e,t),{data:a,reload:t}}const au=`Hello, I am zorrooz, currently focused on structural analysis and functional studies of biological macromolecules, while also exploring computational biology.
`,eu=[{year:"2021 – 2025",title:"Lanzhou University · B.S. in Bioinformatics",desc:"Basic programming training and foundational knowledge in biology"},{year:"2025 – present",title:"Lanzhou University · M.S. in Biology",desc:"Structural and functional studies of biological macromolecules"}],tu=[{title:"Common Tools",items:[{name:"Python · R · Bash · Git",desc:"Primary languages and version control tools for daily research and development"},{name:"Nextflow · Snakemake",desc:"Bioinformatics pipeline orchestration and workflow management"},{name:"RELION · cryoSPARC",desc:"Cryo-EM single-particle reconstruction and 3D structure validation"},{name:"VS Code · Ubuntu · WSL",desc:"Development environment and terminal workflow"}]},{title:"Fields",items:[{name:"Structural Biology",desc:"Cryo-EM"},{name:"Computational Biology",desc:"Binder design and virtual screening"},{name:"Programming",desc:"R · Python · C · Perl · Bash · JavaScript"}]}],lu=[{label:"Email",value:"zorrooz@163.com",link:"mailto:zorrooz@163.com",icon:"fas fa-envelope"},{label:"GitHub",value:"zorrooz",link:"https://github.com/zorrooz",icon:"fab fa-github"}],hb={introduction:au,experience:eu,section:tu,contacts:lu},db=Object.freeze(Object.defineProperty({__proto__:null,contacts:lu,default:hb,experience:eu,introduction:au,section:tu},Symbol.toStringTag,{value:"Module"})),ou=`你好，我是 zorrooz，当前专注于生物大分子的结构解析与功能研究，同时涉猎计算生物学。
`,pu=[{year:"2021 – 2025",title:"兰州大学 本科 · 生物信息学",desc:"基本编程训练与生物学基础知识"},{year:"2025 – 至今",title:"兰州大学 硕士 · 生物学",desc:"生物大分子的结构与功能研究"}],cu=[{title:"常用工具",items:[{name:"Python · R · Bash · Git",desc:"日常研究与开发的主力语言与版本控制工具"},{name:"Nextflow · Snakemake",desc:"生信流水线编排与工作流管理"},{name:"RELION · cryoSPARC",desc:"冷冻电镜单颗粒重构与三维结构验证"},{name:"VS Code · Ubuntu · WSL",desc:"开发环境与终端工作流"}]},{title:"领域",items:[{name:"结构生物学",desc:"冷冻电镜"},{name:"计算生物学",desc:"binder 设计与虚拟筛选"},{name:"编程",desc:"R · Python · C · Perl · Bash · JavaScript"}]}],ru=[{label:"邮箱",value:"zorrooz@163.com",link:"mailto:zorrooz@163.com",icon:"fas fa-envelope"},{label:"GitHub",value:"zorrooz",link:"https://github.com/zorrooz",icon:"fab fa-github"}],mb={introduction:ou,experience:pu,section:cu,contacts:ru},jb=Object.freeze(Object.defineProperty({__proto__:null,contacts:ru,default:mb,experience:pu,introduction:ou,section:cu},Symbol.toStringTag,{value:"Module"})),gb=JSON.parse('[{"title":"Notes","items":[{"name":"ComputationalBiology","title":"Computational Biology","desc":"Mainstream tools for protein design and virtual screening","tags":["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology","Virtual Screening","Molecular Docking","AutoDock Vina"],"stats":{"postsCount":2,"totalWords":1211,"latestDate":"2026-08-04"},"categories":[{"key":"protein-design","title":"Protein Design","articles":[{"title":"Mainstream Tools for Protein Design","articleUrl":"/article/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en","wordCount":618,"tags":["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology"]}],"stats":{"postsCount":1,"totalWords":618,"latestDate":"2026-08-04"}},{"key":"virtual-screening","title":"Virtual Screening","articles":[{"title":"Mainstream Tools for Virtual Screening","articleUrl":"/article/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools-en","wordCount":593,"tags":["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"]}],"stats":{"postsCount":1,"totalWords":593,"latestDate":"2026-08-04"}}],"root":"/article/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en"},{"name":"Programming","title":"Programming","desc":"Detailed tutorials on R, Python, Linux, and Bash programming","tags":["Bash","Shell","Scripting","Tutorial","Linux","Command Line","Python","Advanced","Getting Started","NumPy","Pandas","Data Processing","R Language","ggplot2","Data Visualization","tidyverse","dplyr"],"stats":{"postsCount":8,"totalWords":1972,"latestDate":"2026-08-04"},"categories":[{"key":"python","title":"Python","articles":[{"title":"Python Advanced: Functions, Classes, and Modules","articleUrl":"/article/notes/Programming/python/python-advanced/python-advanced-en","wordCount":244,"tags":["Python","Advanced","Tutorial"]},{"title":"Introduction to Python Programming: Environment, Syntax, and Data Types","articleUrl":"/article/notes/Programming/python/python-basics/python-basics-en","wordCount":333,"tags":["Python","Getting Started","Tutorial"]},{"title":"Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas","articleUrl":"/article/notes/Programming/python/python-data/python-data-en","wordCount":314,"tags":["Python","NumPy","Pandas","Data Processing","Tutorial"]}],"stats":{"postsCount":3,"totalWords":891,"latestDate":"2026-08-04"}},{"key":"r","title":"R","articles":[{"title":"R Language Primer: Data Structures, Vectorization, and Functions","articleUrl":"/article/notes/Programming/r/r-basics/r-basics-en","wordCount":256,"tags":["R Language","Getting Started","Tutorial"]},{"title":"R ggplot2 Data Visualization: Grammar of Graphics and Common Charts","articleUrl":"/article/notes/Programming/r/r-ggplot2/r-ggplot2-en","wordCount":157,"tags":["R Language","ggplot2","Data Visualization","Tutorial"]},{"title":"R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning","articleUrl":"/article/notes/Programming/r/r-tidyverse/r-tidyverse-en","wordCount":245,"tags":["R Language","tidyverse","dplyr","Tutorial"]}],"stats":{"postsCount":3,"totalWords":658,"latestDate":"2026-08-04"}},{"key":"linux","title":"Linux","articles":[{"title":"Linux Command Line Basics: File System, Permissions, and Text Processing","articleUrl":"/article/notes/Programming/linux/linux-basics/linux-basics-en","wordCount":238,"tags":["Linux","Command Line","Tutorial"]}],"stats":{"postsCount":1,"totalWords":238,"latestDate":"2026-08-04"}},{"key":"bash","title":"Bash","articles":[{"title":"Bash Programming: Variables, Conditionals, Loops, and Practical Scripts","articleUrl":"/article/notes/Programming/bash/bash-scripting/bash-scripting-en","wordCount":185,"tags":["Bash","Shell","Scripting","Tutorial"]}],"stats":{"postsCount":1,"totalWords":185,"latestDate":"2026-08-04"}}],"root":"/article/notes/Programming/python/python-advanced/python-advanced-en"},{"name":"StructuralBiology","title":"Structural Biology","desc":"Cryo-EM data processing, visualization of biomacromolecular structures, and atomic modeling","tags":["Cryo-Electron Microscopy","cryo-EM","Review","Data Processing","RELION","Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial","Structure Visualization","PyMOL","ChimeraX"],"stats":{"postsCount":4,"totalWords":2670,"latestDate":"2026-08-04"},"categories":[{"key":"cryoem","title":"Cryo-EM","articles":[{"title":"Review of Cryo-EM Technology","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en","wordCount":977,"tags":["Cryo-Electron Microscopy","cryo-EM","Review"]},{"title":"Cryo-EM Single Particle Analysis: Full Data Processing Workflow","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow-en","wordCount":765,"tags":["Cryo-Electron Microscopy","cryo-EM","Data Processing","RELION"]}],"stats":{"postsCount":2,"totalWords":1742,"latestDate":"2026-08-04"}},{"key":"visualization","title":"Structural Visualization","articles":[{"title":"Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice","articleUrl":"/article/notes/StructuralBiology/visualization/structure-visualization/structure-visualization-en","wordCount":362,"tags":["Structure Visualization","PyMOL","ChimeraX","Tutorial"]}],"stats":{"postsCount":1,"totalWords":362,"latestDate":"2026-08-04"}},{"key":"modeling","title":"Atomic Modeling","articles":[{"title":"Atomic Modeling and Refinement","articleUrl":"/article/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling-en","wordCount":566,"tags":["Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial"]}],"stats":{"postsCount":1,"totalWords":566,"latestDate":"2026-08-04"}}],"root":"/article/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en"}]},{"title":"Projects","items":[{"name":"Plotit","title":"plotit","desc":"Declarative plotting package for R — a verb-prefix API based on ggplot2, in early development stage","tags":["R","ggplot2","Data Visualization","Declarative API"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-12-30"},"categories":[],"root":"","github":"https://github.com/zorrooz/plotit","status":"active","language":"R","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/plotit"},{"name":"ZorroozBlog","title":"zorrooz.github.io","desc":"Personal blog system — built with Vue 3 + Vite, supporting bilingual Chinese/English, theme switching, and Markdown article management","tags":["Vue","Vite","Blog","Markdown"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-08-29"},"categories":[],"root":"","github":"https://github.com/zorrooz/zorrooz.github.io","status":"active","language":"Vue","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/"}]},{"title":"Topics","items":[{"name":"StructureOfProteinDemo","title":"Protein Structure Determination Demo Project","desc":"Protein structure placeholder demo","tags":["Protein","Structure","Cryo-EM"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2026-04-12"},"categories":[],"root":"","doi":"10.1234/demo-structure.2026","status":"completed","journal":"Journal of Demo Science","year":2026,"authors":["Demo Author A","Demo Author B"],"url":"https://www.demo-structure.org"}]}]'),fb=Object.freeze(Object.defineProperty({__proto__:null,default:gb},Symbol.toStringTag,{value:"Module"})),bb=JSON.parse('[{"title":"笔记","items":[{"name":"ComputationalBiology","title":"计算生物学","desc":"蛋白质设计与虚拟筛选的主流工具","tags":["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学","虚拟筛选","分子对接","AutoDock Vina"],"stats":{"postsCount":2,"totalWords":1899,"latestDate":"2026-08-04"},"categories":[{"key":"protein-design","title":"蛋白质设计","articles":[{"title":"蛋白质设计主流工具","articleUrl":"/article/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools","wordCount":1002,"tags":["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学"]}],"stats":{"postsCount":1,"totalWords":1002,"latestDate":"2026-08-04"}},{"key":"virtual-screening","title":"虚拟筛选","articles":[{"title":"虚拟筛选主流工具","articleUrl":"/article/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools","wordCount":897,"tags":["虚拟筛选","分子对接","AutoDock Vina","计算生物学"]}],"stats":{"postsCount":1,"totalWords":897,"latestDate":"2026-08-04"}}],"root":"/article/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools"},{"name":"Programming","title":"编程","desc":"R、Python、Linux 与 Bash 编程的详细教程","tags":["Bash","Shell","脚本编程","教程","Linux","命令行","Python","进阶","入门","NumPy","Pandas","数据处理","R语言","ggplot2","数据可视化","tidyverse","dplyr"],"stats":{"postsCount":8,"totalWords":2903,"latestDate":"2026-08-04"},"categories":[{"key":"python","title":"Python","articles":[{"title":"Python 进阶：函数、类与模块","articleUrl":"/article/notes/Programming/python/python-advanced/python-advanced","wordCount":376,"tags":["Python","进阶","教程"]},{"title":"Python 编程入门：环境、语法与数据类型","articleUrl":"/article/notes/Programming/python/python-basics/python-basics","wordCount":495,"tags":["Python","入门","教程"]},{"title":"Python 数据处理实战：文件 IO、正则与 NumPy/Pandas","articleUrl":"/article/notes/Programming/python/python-data/python-data","wordCount":473,"tags":["Python","NumPy","Pandas","数据处理","教程"]}],"stats":{"postsCount":3,"totalWords":1344,"latestDate":"2026-08-04"}},{"key":"r","title":"R","articles":[{"title":"R 语言入门：数据结构、向量化与函数","articleUrl":"/article/notes/Programming/r/r-basics/r-basics","wordCount":403,"tags":["R语言","入门","教程"]},{"title":"R ggplot2 数据可视化：图层语法与常用图表","articleUrl":"/article/notes/Programming/r/r-ggplot2/r-ggplot2","wordCount":245,"tags":["R语言","ggplot2","数据可视化","教程"]},{"title":"R tidyverse 数据操作：dplyr 管道与数据清洗","articleUrl":"/article/notes/Programming/r/r-tidyverse/r-tidyverse","wordCount":324,"tags":["R语言","tidyverse","dplyr","教程"]}],"stats":{"postsCount":3,"totalWords":972,"latestDate":"2026-08-04"}},{"key":"linux","title":"Linux","articles":[{"title":"Linux 命令行基础：文件系统、权限与文本处理","articleUrl":"/article/notes/Programming/linux/linux-basics/linux-basics","wordCount":328,"tags":["Linux","命令行","教程"]}],"stats":{"postsCount":1,"totalWords":328,"latestDate":"2026-08-04"}},{"key":"bash","title":"Bash","articles":[{"title":"Bash 编程：变量、条件、循环与实用脚本","articleUrl":"/article/notes/Programming/bash/bash-scripting/bash-scripting","wordCount":259,"tags":["Bash","Shell","脚本编程","教程"]}],"stats":{"postsCount":1,"totalWords":259,"latestDate":"2026-08-04"}}],"root":"/article/notes/Programming/python/python-advanced/python-advanced"},{"name":"StructuralBiology","title":"结构生物学","desc":"冷冻电镜数据处理、生物大分子结构可视化与原子建模","tags":["冷冻电镜","cryo-EM","综述","数据处理","RELION","原子建模","Coot","phenix","结构精修","教程","结构可视化","PyMOL","ChimeraX"],"stats":{"postsCount":4,"totalWords":3960,"latestDate":"2026-08-04"},"categories":[{"key":"cryoem","title":"冷冻电镜","articles":[{"title":"冷冻电镜技术综述","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview","wordCount":1485,"tags":["冷冻电镜","cryo-EM","综述"]},{"title":"冷冻电镜单颗粒分析：数据处理全流程","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow","wordCount":1111,"tags":["冷冻电镜","cryo-EM","数据处理","RELION"]}],"stats":{"postsCount":2,"totalWords":2596,"latestDate":"2026-08-04"}},{"key":"visualization","title":"结构可视化","articles":[{"title":"生物大分子结构可视化：PyMOL 与 ChimeraX 实战","articleUrl":"/article/notes/StructuralBiology/visualization/structure-visualization/structure-visualization","wordCount":529,"tags":["结构可视化","PyMOL","ChimeraX","教程"]}],"stats":{"postsCount":1,"totalWords":529,"latestDate":"2026-08-04"}},{"key":"modeling","title":"原子建模","articles":[{"title":"原子建模与精修","articleUrl":"/article/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling","wordCount":835,"tags":["原子建模","Coot","phenix","结构精修","教程"]}],"stats":{"postsCount":1,"totalWords":835,"latestDate":"2026-08-04"}}],"root":"/article/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview"}]},{"title":"项目","items":[{"name":"Plotit","title":"plotit","desc":"R 声明式绘图包——基于 ggplot2 的动词前缀 API，早期开发阶段","tags":["R","ggplot2","数据可视化","声明式API"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-12-30"},"categories":[],"root":"","github":"https://github.com/zorrooz/plotit","status":"active","language":"R","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/plotit"},{"name":"ZorroozBlog","title":"zorrooz.github.io","desc":"个人博客系统——基于 Vue 3 + Vite 构建，支持中英双语、主题切换与 Markdown 文章管理","tags":["Vue","Vite","博客","Markdown"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-08-29"},"categories":[],"root":"","github":"https://github.com/zorrooz/zorrooz.github.io","status":"active","language":"Vue","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/"}]},{"title":"课题","items":[{"name":"StructureOfProteinDemo","title":"蛋白质结构解析演示课题","desc":"蛋白质结构占位demo","tags":["Protein","Structure","Cryo-EM"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2026-04-12"},"categories":[],"root":"","doi":"10.1234/demo-structure.2026","status":"completed","journal":"Journal of Demo Science","year":2026,"authors":["Demo Author A","Demo Author B"],"url":"https://www.demo-structure.org"}]}]'),_b=Object.freeze(Object.defineProperty({__proto__:null,default:bb},Symbol.toStringTag,{value:"Module"})),yb=[{title:"Mainstream Tools for Protein Design",date:"2026-08-04",tags:["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology"],draft:!1,description:"A panorama of protein design tools: Rosetta energy optimization, AlphaFold structure prediction, RFdiffusion generative design, ESMFold fast prediction, and binder design workflow",relativePath:"ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en",wordCount:618,tagCount:5},{title:"Mainstream Tools for Virtual Screening",date:"2026-08-04",tags:["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"],draft:!1,description:"A comprehensive overview of virtual screening tools: docking software such as AutoDock Vina, Glide, and DOCK, receptor/ligand preparation, scoring functions, high-throughput screening workflows, and validation strategies",relativePath:"ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools-en",wordCount:593,tagCount:4},{title:"Bash Programming: Variables, Conditionals, Loops, and Practical Scripts",date:"2026-08-04",tags:["Bash","Shell","Scripting","Tutorial"],draft:!1,description:"A complete Bash scripting tutorial: variables and parameters, conditional statements, loops, functions, arrays, and string processing, with practical scripts for batch processing and bioinformatics",relativePath:"Programming/bash/bash-scripting/bash-scripting-en",wordCount:185,tagCount:4},{title:"Linux Command Line Basics: File System, Permissions, and Text Processing",date:"2026-08-04",tags:["Linux","Command Line","Tutorial"],draft:!1,description:"Core skills for the Linux command line: file system navigation, file operations, permission management, the three musketeers of text processing (grep/sed/awk), and process management",relativePath:"Programming/linux/linux-basics/linux-basics-en",wordCount:238,tagCount:3},{title:"Python Advanced: Functions, Classes, and Modules",date:"2026-08-04",tags:["Python","Advanced","Tutorial"],draft:!1,description:"Deep dive into Python functional and object-oriented programming: parameter passing, lambda, decorators, classes and inheritance, exception handling, modules and packages",relativePath:"Programming/python/python-advanced/python-advanced-en",wordCount:244,tagCount:3},{title:"Introduction to Python Programming: Environment, Syntax, and Data Types",date:"2026-08-04",tags:["Python","Getting Started","Tutorial"],draft:!1,description:"Learn Python from scratch: environment setup, basic syntax, built-in data types, and control flow, with runnable example code.",relativePath:"Programming/python/python-basics/python-basics-en",wordCount:333,tagCount:3},{title:"Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas",date:"2026-08-04",tags:["Python","NumPy","Pandas","Data Processing","Tutorial"],draft:!1,description:"A complete practical tutorial on file I/O, regular expressions, and the scientific computing stack (NumPy/Pandas), covering common data scenarios in bioinformatics",relativePath:"Programming/python/python-data/python-data-en",wordCount:314,tagCount:5},{title:"R Language Primer: Data Structures, Vectorization, and Functions",date:"2026-08-04",tags:["R Language","Getting Started","Tutorial"],draft:!1,description:"Core R fundamentals: environment setup, vectors/matrices/data frames/lists, vectorized operations, control flow, and function definitions",relativePath:"Programming/r/r-basics/r-basics-en",wordCount:256,tagCount:3},{title:"R ggplot2 Data Visualization: Grammar of Graphics and Common Charts",date:"2026-08-04",tags:["R Language","ggplot2","Data Visualization","Tutorial"],draft:!1,description:"Introduction to the ggplot2 Grammar of Graphics: scatter plots, boxplots, histograms, bar charts, and theme customization for publication-ready figures",relativePath:"Programming/r/r-ggplot2/r-ggplot2-en",wordCount:157,tagCount:4},{title:"R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning",date:"2026-08-04",tags:["R Language","tidyverse","dplyr","Tutorial"],draft:!1,description:"Chain dplyr verbs (filter/select/mutate/group_by/summarise/join) with pipes, and use tidyr to complete a full data cleaning workflow.",relativePath:"Programming/r/r-tidyverse/r-tidyverse-en",wordCount:245,tagCount:4},{title:"Review of Cryo-EM Technology",date:"2026-08-04",tags:["Cryo-Electron Microscopy","cryo-EM","Review"],draft:!1,description:"The technological revolution of single-particle cryo-EM: direct electron detectors, the stability revolution, and the advent of the AI era, along with the future landscape of structural biology",relativePath:"StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en",wordCount:977,tagCount:3},{title:"Cryo-EM Single Particle Analysis: Full Data Processing Workflow",date:"2026-08-04",tags:["Cryo-Electron Microscopy","cryo-EM","Data Processing","RELION"],draft:!1,description:"From micrographs to atomic models: the complete workflow of cryo-EM single particle analysis (SPA), covering preprocessing, particle picking, 2D/3D classification, refinement, and resolution assessment",relativePath:"StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow-en",wordCount:765,tagCount:4},{title:"Atomic Modeling and Refinement",date:"2026-08-04",tags:["Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial"],draft:!1,description:"Complete workflow for atomic model building and refinement from cryo-EM density maps: automatic modeling, Coot manual correction, phenix refinement, model quality assessment metrics",relativePath:"StructuralBiology/modeling/atomic-modeling/atomic-modeling-en",wordCount:566,tagCount:5},{title:"Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice",date:"2026-08-04",tags:["Structure Visualization","PyMOL","ChimeraX","Tutorial"],draft:!1,description:"Visualizing protein/nucleic acid structures with PyMOL and UCSF ChimeraX: PDB data retrieval, rendering modes, script batch processing, and cryo-EM density map display",relativePath:"StructuralBiology/visualization/structure-visualization/structure-visualization-en",wordCount:362,tagCount:4}],vb=Object.freeze(Object.defineProperty({__proto__:null,default:yb},Symbol.toStringTag,{value:"Module"})),wb=[{title:"蛋白质设计主流工具",date:"2026-08-04",tags:["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学"],draft:!1,description:"蛋白质设计工具全景：Rosetta 能量优化、AlphaFold 结构预测、RFdiffusion 生成式设计、ESMFold 快速预测与 binder 设计工作流",relativePath:"ComputationalBiology/protein-design/protein-design-tools/protein-design-tools",wordCount:1002,tagCount:5},{title:"虚拟筛选主流工具",date:"2026-08-04",tags:["虚拟筛选","分子对接","AutoDock Vina","计算生物学"],draft:!1,description:"虚拟筛选工具全景：AutoDock Vina、Glide、DOCK 等对接软件，受体/配体准备、打分函数、高通量筛选流程与验证策略",relativePath:"ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools",wordCount:897,tagCount:4},{title:"Bash 编程：变量、条件、循环与实用脚本",date:"2026-08-04",tags:["Bash","Shell","脚本编程","教程"],draft:!1,description:"Bash 脚本编程完整教程：变量与参数、条件判断、循环、函数、数组与字符串处理，附批处理与生信实用脚本",relativePath:"Programming/bash/bash-scripting/bash-scripting",wordCount:259,tagCount:4},{title:"Linux 命令行基础：文件系统、权限与文本处理",date:"2026-08-04",tags:["Linux","命令行","教程"],draft:!1,description:"Linux 命令行核心技能：文件系统导航、文件操作、权限管理、文本处理三剑客（grep/sed/awk）与进程管理",relativePath:"Programming/linux/linux-basics/linux-basics",wordCount:328,tagCount:3},{title:"Python 进阶：函数、类与模块",date:"2026-08-04",tags:["Python","进阶","教程"],draft:!1,description:"深入 Python 函数式与面向对象编程：参数传递、lambda、装饰器、类与继承、异常处理、模块与包",relativePath:"Programming/python/python-advanced/python-advanced",wordCount:376,tagCount:3},{title:"Python 编程入门：环境、语法与数据类型",date:"2026-08-04",tags:["Python","入门","教程"],draft:!1,description:"从零开始学习 Python：环境搭建、基础语法、内置数据类型与流程控制，附带可运行的示例代码",relativePath:"Programming/python/python-basics/python-basics",wordCount:495,tagCount:3},{title:"Python 数据处理实战：文件 IO、正则与 NumPy/Pandas",date:"2026-08-04",tags:["Python","NumPy","Pandas","数据处理","教程"],draft:!1,description:"文件读写、正则表达式与科学计算三件套（NumPy/Pandas）的完整实战教程，覆盖生物信息学常见数据场景",relativePath:"Programming/python/python-data/python-data",wordCount:473,tagCount:5},{title:"R 语言入门：数据结构、向量化与函数",date:"2026-08-04",tags:["R语言","入门","教程"],draft:!1,description:"R 语言核心基础：环境配置、向量/矩阵/数据框/列表、向量化运算、控制流与函数定义",relativePath:"Programming/r/r-basics/r-basics",wordCount:403,tagCount:3},{title:"R ggplot2 数据可视化：图层语法与常用图表",date:"2026-08-04",tags:["R语言","ggplot2","数据可视化","教程"],draft:!1,description:"ggplot2 图层语法入门：散点图、箱线图、直方图、柱状图与主题美化，输出出版级图表",relativePath:"Programming/r/r-ggplot2/r-ggplot2",wordCount:245,tagCount:4},{title:"R tidyverse 数据操作：dplyr 管道与数据清洗",date:"2026-08-04",tags:["R语言","tidyverse","dplyr","教程"],draft:!1,description:"用管道符串联 dplyr 动词（filter/select/mutate/group_by/summarise/join），配合 tidyr 完成完整的数据清洗流程",relativePath:"Programming/r/r-tidyverse/r-tidyverse",wordCount:324,tagCount:4},{title:"冷冻电镜技术综述",date:"2026-08-04",tags:["冷冻电镜","cryo-EM","综述"],draft:!1,description:"冷冻电镜单颗粒分析的技术革命：直接电子探测器、稳定性革命与 AI 时代的到来，以及结构生物学的未来图景",relativePath:"StructuralBiology/cryoem/cryoem-overview/cryoem-overview",wordCount:1485,tagCount:3},{title:"冷冻电镜单颗粒分析：数据处理全流程",date:"2026-08-04",tags:["冷冻电镜","cryo-EM","数据处理","RELION"],draft:!1,description:"从显微照片到原子模型：冷冻电镜单颗粒分析（SPA）完整流程，涵盖预处理、颗粒挑选、2D/3D 分类、精修与分辨率评估",relativePath:"StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow",wordCount:1111,tagCount:4},{title:"原子建模与精修",date:"2026-08-04",tags:["原子建模","Coot","phenix","结构精修","教程"],draft:!1,description:"冷冻电镜密度图的原子模型搭建与精修全流程：自动建模、Coot 手动修正、phenix 精修、模型质量评估指标",relativePath:"StructuralBiology/modeling/atomic-modeling/atomic-modeling",wordCount:835,tagCount:5},{title:"生物大分子结构可视化：PyMOL 与 ChimeraX 实战",date:"2026-08-04",tags:["结构可视化","PyMOL","ChimeraX","教程"],draft:!1,description:"用 PyMOL 与 UCSF ChimeraX 可视化蛋白质/核酸结构：PDB 数据获取、渲染模式、脚本批处理与冷冻电镜密度图展示",relativePath:"StructuralBiology/visualization/structure-visualization/structure-visualization",wordCount:529,tagCount:4}],Cb=Object.freeze(Object.defineProperty({__proto__:null,default:wb},Symbol.toStringTag,{value:"Module"})),kb=[{id:1,no:1,title:"Mainstream Tools for Protein Design",date:"2026-08-04",category:["ComputationalBiology","protein-design"],tags:["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology"],preview:"A panorama of protein design tools: Rosetta energy optimization, AlphaFold structure prediction, RFdiffusion generative design, ESMFold fast prediction, and binder design workflow",wordCount:618},{id:2,no:2,title:"Mainstream Tools for Virtual Screening",date:"2026-08-04",category:["ComputationalBiology","virtual-screening"],tags:["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"],preview:"A comprehensive overview of virtual screening tools: docking software such as AutoDock Vina, Glide, and DOCK, receptor/ligand preparation, scoring functions, high-throughput screening workflows, and validation strategies",wordCount:593},{id:3,no:3,title:"Bash Programming: Variables, Conditionals, Loops, and Practical Scripts",date:"2026-08-04",category:["Programming","bash"],tags:["Bash","Shell","Scripting","Tutorial"],preview:"A complete Bash scripting tutorial: variables and parameters, conditional statements, loops, functions, arrays, and string processing, with practical scripts for batch processing and bioinformatics",wordCount:185},{id:4,no:4,title:"Linux Command Line Basics: File System, Permissions, and Text Processing",date:"2026-08-04",category:["Programming","linux"],tags:["Linux","Command Line","Tutorial"],preview:"Core skills for the Linux command line: file system navigation, file operations, permission management, the three musketeers of text processing (grep/sed/awk), and process management",wordCount:238},{id:5,no:5,title:"Python Advanced: Functions, Classes, and Modules",date:"2026-08-04",category:["Programming","python"],tags:["Python","Advanced","Tutorial"],preview:"Deep dive into Python functional and object-oriented programming: parameter passing, lambda, decorators, classes and inheritance, exception handling, modules and packages",wordCount:244},{id:6,no:6,title:"Introduction to Python Programming: Environment, Syntax, and Data Types",date:"2026-08-04",category:["Programming","python"],tags:["Python","Getting Started","Tutorial"],preview:"Learn Python from scratch: environment setup, basic syntax, built-in data types, and control flow, with runnable example code.",wordCount:333},{id:7,no:7,title:"Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas",date:"2026-08-04",category:["Programming","python"],tags:["Python","NumPy","Pandas","Data Processing","Tutorial"],preview:"A complete practical tutorial on file I/O, regular expressions, and the scientific computing stack (NumPy/Pandas), covering common data scenarios in bioinformatics",wordCount:314},{id:8,no:8,title:"R Language Primer: Data Structures, Vectorization, and Functions",date:"2026-08-04",category:["Programming","r"],tags:["R Language","Getting Started","Tutorial"],preview:"Core R fundamentals: environment setup, vectors/matrices/data frames/lists, vectorized operations, control flow, and function definitions",wordCount:256},{id:9,no:9,title:"R ggplot2 Data Visualization: Grammar of Graphics and Common Charts",date:"2026-08-04",category:["Programming","r"],tags:["R Language","ggplot2","Data Visualization","Tutorial"],preview:"Introduction to the ggplot2 Grammar of Graphics: scatter plots, boxplots, histograms, bar charts, and theme customization for publication-ready figures",wordCount:157},{id:10,no:10,title:"R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning",date:"2026-08-04",category:["Programming","r"],tags:["R Language","tidyverse","dplyr","Tutorial"],preview:"Chain dplyr verbs (filter/select/mutate/group_by/summarise/join) with pipes, and use tidyr to complete a full data cleaning workflow.",wordCount:245},{id:11,no:11,title:"Review of Cryo-EM Technology",date:"2026-08-04",category:["StructuralBiology","cryoem"],tags:["Cryo-Electron Microscopy","cryo-EM","Review"],preview:"The technological revolution of single-particle cryo-EM: direct electron detectors, the stability revolution, and the advent of the AI era, along with the future landscape of structural biology",wordCount:977},{id:12,no:12,title:"Cryo-EM Single Particle Analysis: Full Data Processing Workflow",date:"2026-08-04",category:["StructuralBiology","cryoem"],tags:["Cryo-Electron Microscopy","cryo-EM","Data Processing","RELION"],preview:"From micrographs to atomic models: the complete workflow of cryo-EM single particle analysis (SPA), covering preprocessing, particle picking, 2D/3D classification, refinement, and resolution assessment",wordCount:765},{id:13,no:13,title:"Atomic Modeling and Refinement",date:"2026-08-04",category:["StructuralBiology","modeling"],tags:["Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial"],preview:"Complete workflow for atomic model building and refinement from cryo-EM density maps: automatic modeling, Coot manual correction, phenix refinement, model quality assessment metrics",wordCount:566},{id:14,no:14,title:"Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice",date:"2026-08-04",category:["StructuralBiology","visualization"],tags:["Structure Visualization","PyMOL","ChimeraX","Tutorial"],preview:"Visualizing protein/nucleic acid structures with PyMOL and UCSF ChimeraX: PDB data retrieval, rendering modes, script batch processing, and cryo-EM density map display",wordCount:362}],xb=Object.freeze(Object.defineProperty({__proto__:null,default:kb},Symbol.toStringTag,{value:"Module"})),Sb=[{id:1,no:1,title:"蛋白质设计主流工具",date:"2026-08-04",category:["计算生物学","protein-design"],tags:["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学"],preview:"蛋白质设计工具全景：Rosetta 能量优化、AlphaFold 结构预测、RFdiffusion 生成式设计、ESMFold 快速预测与 binder 设计工作流",wordCount:1002},{id:2,no:2,title:"虚拟筛选主流工具",date:"2026-08-04",category:["计算生物学","virtual-screening"],tags:["虚拟筛选","分子对接","AutoDock Vina","计算生物学"],preview:"虚拟筛选工具全景：AutoDock Vina、Glide、DOCK 等对接软件，受体/配体准备、打分函数、高通量筛选流程与验证策略",wordCount:897},{id:3,no:3,title:"Bash 编程：变量、条件、循环与实用脚本",date:"2026-08-04",category:["编程","bash"],tags:["Bash","Shell","脚本编程","教程"],preview:"Bash 脚本编程完整教程：变量与参数、条件判断、循环、函数、数组与字符串处理，附批处理与生信实用脚本",wordCount:259},{id:4,no:4,title:"Linux 命令行基础：文件系统、权限与文本处理",date:"2026-08-04",category:["编程","linux"],tags:["Linux","命令行","教程"],preview:"Linux 命令行核心技能：文件系统导航、文件操作、权限管理、文本处理三剑客（grep/sed/awk）与进程管理",wordCount:328},{id:5,no:5,title:"Python 进阶：函数、类与模块",date:"2026-08-04",category:["编程","python"],tags:["Python","进阶","教程"],preview:"深入 Python 函数式与面向对象编程：参数传递、lambda、装饰器、类与继承、异常处理、模块与包",wordCount:376},{id:6,no:6,title:"Python 编程入门：环境、语法与数据类型",date:"2026-08-04",category:["编程","python"],tags:["Python","入门","教程"],preview:"从零开始学习 Python：环境搭建、基础语法、内置数据类型与流程控制，附带可运行的示例代码",wordCount:495},{id:7,no:7,title:"Python 数据处理实战：文件 IO、正则与 NumPy/Pandas",date:"2026-08-04",category:["编程","python"],tags:["Python","NumPy","Pandas","数据处理","教程"],preview:"文件读写、正则表达式与科学计算三件套（NumPy/Pandas）的完整实战教程，覆盖生物信息学常见数据场景",wordCount:473},{id:8,no:8,title:"R 语言入门：数据结构、向量化与函数",date:"2026-08-04",category:["编程","r"],tags:["R语言","入门","教程"],preview:"R 语言核心基础：环境配置、向量/矩阵/数据框/列表、向量化运算、控制流与函数定义",wordCount:403},{id:9,no:9,title:"R ggplot2 数据可视化：图层语法与常用图表",date:"2026-08-04",category:["编程","r"],tags:["R语言","ggplot2","数据可视化","教程"],preview:"ggplot2 图层语法入门：散点图、箱线图、直方图、柱状图与主题美化，输出出版级图表",wordCount:245},{id:10,no:10,title:"R tidyverse 数据操作：dplyr 管道与数据清洗",date:"2026-08-04",category:["编程","r"],tags:["R语言","tidyverse","dplyr","教程"],preview:"用管道符串联 dplyr 动词（filter/select/mutate/group_by/summarise/join），配合 tidyr 完成完整的数据清洗流程",wordCount:324},{id:11,no:11,title:"冷冻电镜技术综述",date:"2026-08-04",category:["结构生物学","cryoem"],tags:["冷冻电镜","cryo-EM","综述"],preview:"冷冻电镜单颗粒分析的技术革命：直接电子探测器、稳定性革命与 AI 时代的到来，以及结构生物学的未来图景",wordCount:1485},{id:12,no:12,title:"冷冻电镜单颗粒分析：数据处理全流程",date:"2026-08-04",category:["结构生物学","cryoem"],tags:["冷冻电镜","cryo-EM","数据处理","RELION"],preview:"从显微照片到原子模型：冷冻电镜单颗粒分析（SPA）完整流程，涵盖预处理、颗粒挑选、2D/3D 分类、精修与分辨率评估",wordCount:1111},{id:13,no:13,title:"原子建模与精修",date:"2026-08-04",category:["结构生物学","modeling"],tags:["原子建模","Coot","phenix","结构精修","教程"],preview:"冷冻电镜密度图的原子模型搭建与精修全流程：自动建模、Coot 手动修正、phenix 精修、模型质量评估指标",wordCount:835},{id:14,no:14,title:"生物大分子结构可视化：PyMOL 与 ChimeraX 实战",date:"2026-08-04",category:["结构生物学","visualization"],tags:["结构可视化","PyMOL","ChimeraX","教程"],preview:"用 PyMOL 与 UCSF ChimeraX 可视化蛋白质/核酸结构：PDB 数据获取、渲染模式、脚本批处理与冷冻电镜密度图展示",wordCount:529}],Pb=Object.freeze(Object.defineProperty({__proto__:null,default:Sb},Symbol.toStringTag,{value:"Module"})),Eb=[{name:"Plotit",title:"plotit",desc:"Declarative plotting package for R — a verb-prefix API based on ggplot2, in early development stage",date:"2025-12-30",tags:["R","ggplot2","Data Visualization","Declarative API"],github:"https://github.com/zorrooz/plotit",url:"https://zorrooz.github.io/plotit",status:"active",language:"R",stars:0,license:"MIT",version:""},{name:"ZorroozBlog",title:"zorrooz.github.io",desc:"Personal blog system — built with Vue 3 + Vite, supporting bilingual Chinese/English, theme switching, and Markdown article management",date:"2025-08-29",tags:["Vue","Vite","Blog","Markdown"],github:"https://github.com/zorrooz/zorrooz.github.io",url:"https://zorrooz.github.io/",status:"active",language:"Vue",stars:0,license:"MIT",version:""}],Tb=Object.freeze(Object.defineProperty({__proto__:null,default:Eb},Symbol.toStringTag,{value:"Module"})),Ab=[{name:"Plotit",title:"plotit",desc:"R 声明式绘图包——基于 ggplot2 的动词前缀 API，早期开发阶段",date:"2025-12-30",tags:["R","ggplot2","数据可视化","声明式API"],github:"https://github.com/zorrooz/plotit",url:"https://zorrooz.github.io/plotit",status:"active",language:"R",stars:0,license:"MIT",version:""},{name:"ZorroozBlog",title:"zorrooz.github.io",desc:"个人博客系统——基于 Vue 3 + Vite 构建，支持中英双语、主题切换与 Markdown 文章管理",date:"2025-08-29",tags:["Vue","Vite","博客","Markdown"],github:"https://github.com/zorrooz/zorrooz.github.io",url:"https://zorrooz.github.io/",status:"active",language:"Vue",stars:0,license:"MIT",version:""}],Rb=Object.freeze(Object.defineProperty({__proto__:null,default:Ab},Symbol.toStringTag,{value:"Module"})),Lb=JSON.parse(`[{"title":"Data Analysis","children":[{"title":"Data Exploration","children":[{"title":"R Language","items":[{"name":"tidyverse","url":"https://www.tidyverse.org/","desc":"A suite of R packages for data science, covering data import, cleaning, transformation, and visualization"},{"name":"data.table","url":"https://rdatatable.gitlab.io/data.table/","desc":"High-performance data manipulation package, extremely fast for processing data frames with millions of rows"},{"name":"readr","url":"https://readr.tidyverse.org/","desc":"A tool for fast reading of CSV/TSV and other tabular files, with automatic column type inference"}]},{"title":"Python","items":[{"name":"NumPy","url":"https://numpy.org/","desc":"The fundamental library for scientific computing in Python, providing multidimensional arrays and vectorized operations"},{"name":"Pandas","url":"https://pandas.pydata.org/","desc":"A powerful data analysis library based on Python, providing efficient data manipulation and processing capabilities"},{"name":"Polars","url":"https://www.pola.rs/","desc":"A high-performance data processing library based on Apache Arrow, supporting multiple language interfaces"},{"name":"Jupyter","url":"https://jupyter.org/","desc":"An interactive notebook environment that supports mixing code, text, and visualizations"}]}]},{"title":"Data Visualization","children":[{"title":"R Language","items":[{"name":"ggplot2","url":"https://ggplot2.tidyverse.org/","desc":"A powerful data visualization package in R, based on the grammar of graphics"},{"name":"plotly (R)","url":"https://plotly.com/r/","desc":"The R interface to the interactive chart library, capable of generating web-interactive graphics"},{"name":"patchwork","url":"https://patchwork.data-imaginist.com/","desc":"A multi-panel composition tool that easily combines ggplot graphics using + and / symbols"},{"name":"RColorBrewer","url":"https://cran.r-project.org/web/packages/RColorBrewer/","desc":"Publication-quality color schemes, providing carefully designed discrete palettes"}]},{"title":"Python","items":[{"name":"Matplotlib","url":"https://matplotlib.org/","desc":"A widely used plotting library in Python, supporting many chart types"},{"name":"Seaborn","url":"https://seaborn.pydata.org/","desc":"A statistical visualization library based on Matplotlib, with built-in attractive themes and statistical plots"},{"name":"Plotly","url":"https://plotly.com/python/","desc":"An interactive visualization library supporting scatter, 3D, maps, and many other chart types"},{"name":"Altair","url":"https://altair-viz.github.io/","desc":"A declarative statistical visualization library based on Vega-Lite"}]}]},{"title":"Statistical Analysis","children":[{"title":"R Language","items":[{"name":"R Stats","url":"https://www.r-project.org/","desc":"The built-in statistical analysis functionality in R, covering a wide range of statistical methods and models"},{"name":"rstatix","url":"https://rpkgs.datanovia.com/rstatix/","desc":"Pipe-friendly wrappers for statistical tests (t-tests, ANOVA, rank-sum tests, etc.)"},{"name":"broom","url":"https://broom.tidymodels.org/","desc":"A tool for tidying statistical model outputs into clean data frames"},{"name":"lme4","url":"https://cran.r-project.org/web/packages/lme4/","desc":"Linear and nonlinear mixed-effects models, suitable for repeated-measures data"}]},{"title":"Python","items":[{"name":"SciPy Stats","url":"https://scipy.org/","desc":"The statistics module in the Python scientific computing library SciPy, providing many statistical distributions and test methods"},{"name":"statsmodels","url":"https://www.statsmodels.org/","desc":"A statistical modeling library supporting regression, time series, hypothesis testing, and more"},{"name":"pingouin","url":"https://pingouin-stats.org/","desc":"A simple and user-friendly statistical testing library covering common parametric and nonparametric tests"}]}]},{"title":"Machine Learning","children":[{"title":"Python","items":[{"name":"scikit-learn","url":"https://scikit-learn.org/","desc":"A machine learning library based on Python, providing a rich set of algorithms and tools"},{"name":"PyTorch","url":"https://pytorch.org/","desc":"A deep learning framework supporting dynamic computation graphs and efficient tensor operations"},{"name":"TensorFlow","url":"https://www.tensorflow.org/","desc":"A deep learning framework with a mature ecosystem, supporting production-level deployment"},{"name":"XGBoost","url":"https://xgboost.readthedocs.io/","desc":"A gradient boosting tree library, the top choice for tabular data competitions and engineering"},{"name":"LightGBM","url":"https://lightgbm.readthedocs.io/","desc":"A gradient boosting framework from Microsoft, with fast training speed and low memory usage"}]},{"title":"R Language","items":[{"name":"tidymodels","url":"https://www.tidymodels.org/","desc":"A unified modeling framework in R, covering data preprocessing, modeling, and evaluation"}]}]}]},{"title":"Omics","children":[{"title":"Data Platforms","children":[{"title":"Sequence Databases","items":[{"name":"NCBI","url":"https://www.ncbi.nlm.nih.gov/","desc":"The National Center for Biotechnology Information in the United States, providing important databases such as GenBank"},{"name":"EMBL-EBI","url":"https://www.ebi.ac.uk/","desc":"The European Molecular Biology Laboratory, providing a variety of bioinformatics resources and tools"},{"name":"DDBJ","url":"https://www.ddbj.nig.ac.jp/","desc":"The DNA Data Bank of Japan, providing storage and access for nucleic acid and protein sequence data"},{"name":"UniProt","url":"https://www.uniprot.org/","desc":"The most comprehensive database for protein sequences and functional annotation"},{"name":"Ensembl","url":"https://www.ensembl.org/","desc":"A genome annotation and comparative genomics database for vertebrates"},{"name":"UCSC Genome Browser","url":"https://genome.ucsc.edu/","desc":"A genome browser supporting visualization of multiple species genomes and custom tracks"}]},{"title":"Sequencing Data","items":[{"name":"SRA","url":"https://www.ncbi.nlm.nih.gov/sra","desc":"The NCBI Sequence Read Archive, storing raw high-throughput sequencing data"},{"name":"GEO","url":"https://www.ncbi.nlm.nih.gov/geo/","desc":"A gene expression database containing microarray and sequencing expression data"},{"name":"ENA","url":"https://www.ebi.ac.uk/ena","desc":"The European Nucleotide Archive, Europe's storage center for raw sequencing data"}]},{"title":"Protein Interactions and Pathways","items":[{"name":"STRING","url":"https://string-db.org/","desc":"A protein-protein interaction database supporting multiple species"},{"name":"KEGG","url":"https://www.kegg.jp/","desc":"A comprehensive database of pathways, genes, and compounds"},{"name":"GO","url":"http://geneontology.org/","desc":"Gene Ontology, providing a standardized system for gene function annotation"}]}]},{"title":"Analysis Software","children":[{"title":"Command Line","items":[{"name":"FastQC","url":"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/","desc":"A quality assessment tool for sequencing data, generating visual QC reports"},{"name":"Trimmomatic","url":"http://www.usadellab.org/cms/?page=trimmomatic","desc":"A read trimming tool for sequencing data, removing adapters and low-quality bases"},{"name":"cutadapt","url":"https://cutadapt.readthedocs.io/","desc":"A command-line adapter removal tool supporting multiple sequencing platforms"},{"name":"BWA","url":"https://github.com/lh3/bwa","desc":"A short-read alignment tool, the standard choice for genome alignment"},{"name":"STAR","url":"https://github.com/alexdobin/STAR","desc":"An ultra-fast splice-aware aligner for RNA-seq"},{"name":"HISAT2","url":"https://daehwankimlab.github.io/hisat2/","desc":"A graph-based index alignment tool for RNA-seq"},{"name":"SAMtools","url":"https://www.htslib.org/","desc":"A toolset for processing SAM/BAM files, the standard for post-alignment processing"},{"name":"bcftools","url":"https://samtools.github.io/bcftools/","desc":"A tool for processing VCF/BCF variant files, supporting filtering and statistics"},{"name":"GATK","url":"https://gatk.broadinstitute.org/","desc":"The Genome Analysis Toolkit, the industry standard for variant detection"},{"name":"MultiQC","url":"https://multiqc.info/","desc":"Aggregates multiple QC reports into a single interactive HTML report"}]},{"title":"Python","items":[{"name":"Biopython","url":"https://biopython.org/","desc":"The core Python library for bioinformatics, covering sequences, structures, and file parsing"},{"name":"scikit-bio","url":"https://scikit-bio.org/","desc":"A Python library for bioinformatics, providing various tools for biological data processing and analysis"},{"name":"scanpy","url":"https://scanpy.readthedocs.io/","desc":"A Python library for single-cell transcriptomics analysis, deeply integrated with AnnData"}]},{"title":"R Language","items":[{"name":"Bioconductor","url":"https://bioconductor.org/","desc":"A collection of bioinformatics software packages based on R, covering genomics, transcriptomics, and other fields"},{"name":"Seurat","url":"https://satijalab.org/seurat/","desc":"A mainstream R package for single-cell transcriptomics analysis, providing a complete analysis workflow"},{"name":"DESeq2","url":"https://bioconductor.org/packages/release/bioc/html/DESeq2.html","desc":"The standard tool for differential expression analysis based on the negative binomial model"},{"name":"edgeR","url":"https://bioconductor.org/packages/release/bioc/html/edgeR.html","desc":"A differential expression analysis tool based on empirical Bayes"},{"name":"limma","url":"https://bioconductor.org/packages/release/bioc/html/limma.html","desc":"A differential analysis tool using linear models, applicable to both microarray and sequencing data"}]}]}]},{"title":"Cryo-EM Structure Determination","children":[{"title":"Software Tools","children":[{"title":"3D Reconstruction","items":[{"name":"RELION","url":"https://www3.mrc-lmb.cam.ac.uk/groups/scheres/index.php","desc":"A Bayesian approach-based software for cryo-EM 3D reconstruction, widely used for high-resolution structure determination"},{"name":"cryoSPARC","url":"https://cryoem.slac.stanford.edu/","desc":"A high-performance cryo-EM data processing platform with an intuitive user interface"},{"name":"cisTEM","url":"https://cistem.org/","desc":"A free, open-source cryo-EM processing software supporting a complete single-particle workflow"},{"name":"SPHIRE","url":"https://sphire.mpg.de/","desc":"A modular cryo-EM processing suite supporting high-resolution reconstruction"}]},{"title":"Preprocessing","items":[{"name":"MotionCor2","url":"https://emcore.ucsf.edu/ucsf-software","desc":"A frame-by-frame motion correction tool that eliminates sample drift caused by the electron beam"},{"name":"CTFFIND4","url":"https://grigoriefflab.umassmed.edu/ctffind4","desc":"A classic tool for CTF estimation, evaluating defocus and astigmatism"},{"name":"Gctf","url":"https://www.mrc-lmb.cam.ac.uk/kzhang/Gctf/","desc":"A GPU-accelerated CTF estimation tool with high speed and precision"},{"name":"Warp","url":"http://www.warpem.com/","desc":"A fast preprocessing tool integrating motion correction, CTF estimation, and particle picking"}]},{"title":"Particle Picking","items":[{"name":"crYOLO","url":"https://cryolo.readthedocs.io/","desc":"A deep learning-based particle picking tool supporting real-time picking"},{"name":"Topaz","url":"https://topaz.readthedocs.io/","desc":"A CNN-based particle picking tool with excellent performance on low signal-to-noise ratio samples"}]},{"title":"Modeling and Visualization","items":[{"name":"PyMOL","url":"https://pymol.org/","desc":"A molecular visualization software, the standard tool for scripted rendering and structure analysis"},{"name":"UCSF ChimeraX","url":"https://www.cgl.ucsf.edu/chimerax/","desc":"A modern molecular visualization software with outstanding density map processing capabilities"},{"name":"Coot","url":"https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/","desc":"A manual model building tool for interactive adjustment of density maps and models"},{"name":"phenix","url":"https://phenix-online.org/","desc":"A structure refinement and validation suite supporting cryo-EM and crystallography"},{"name":"ISOLDE","url":"https://isolde.cimr.cam.ac.uk/","desc":"An interactive real-time molecular dynamics modeling tool based on ChimeraX"},{"name":"ModelAngelo","url":"https://github.com/3dem/model-angelo","desc":"A deep learning-based automatic atomic modeling tool for density maps"}]}]},{"title":"Database Resources","children":[{"title":"Structures and Density Maps","items":[{"name":"EMDB","url":"https://www.ebi.ac.uk/emdb/","desc":"The Electron Microscopy Data Bank, storing and distributing 3D density maps from cryo-EM"},{"name":"EMPIAR","url":"https://www.ebi.ac.uk/empiar/","desc":"The Electron Microscopy Public Image Archive, providing raw movies and particle data"},{"name":"RCSB PDB","url":"https://www.rcsb.org/","desc":"The Protein Data Bank, containing biological macromolecular structures determined by X-ray, NMR, and cryo-EM"},{"name":"wwPDB","url":"https://www.wwpdb.org/","desc":"The worldwide PDB coordination body, uniformly managing structural data standards"},{"name":"PDBe","url":"https://www.ebi.ac.uk/pdbe/","desc":"The European PDB node, providing structural data and visualization tools"},{"name":"AlphaFold DB","url":"https://alphafold.ebi.ac.uk/","desc":"The AlphaFold protein structure database, covering hundreds of millions of proteins"}]}]}]}]`),Mb=Object.freeze(Object.defineProperty({__proto__:null,default:Lb},Symbol.toStringTag,{value:"Module"})),Db=JSON.parse('[{"title":"数据分析","children":[{"title":"数据探索","children":[{"title":"R 语言","items":[{"name":"tidyverse","url":"https://www.tidyverse.org/","desc":"R 语言中用于数据科学的一整套包，涵盖数据导入、清洗、转换和可视化"},{"name":"data.table","url":"https://rdatatable.gitlab.io/data.table/","desc":"高性能数据操作包，百万行级数据框处理速度极快"},{"name":"readr","url":"https://readr.tidyverse.org/","desc":"快速读取 CSV/TSV 等表格文件的工具，自动推断列类型"}]},{"title":"Python","items":[{"name":"NumPy","url":"https://numpy.org/","desc":"Python 科学计算基础库，提供多维数组与向量化运算"},{"name":"Pandas","url":"https://pandas.pydata.org/","desc":"基于 Python 的强大数据分析库，提供高效的数据操作和处理功能"},{"name":"Polars","url":"https://www.pola.rs/","desc":"基于 Apache Arrow 的高性能数据处理库，支持多语言接口"},{"name":"Jupyter","url":"https://jupyter.org/","desc":"交互式笔记本环境，支持代码、文本与可视化混排"}]}]},{"title":"数据可视化","children":[{"title":"R 语言","items":[{"name":"ggplot2","url":"https://ggplot2.tidyverse.org/","desc":"R 语言中功能强大的数据可视化包，基于语法图形学"},{"name":"plotly（R）","url":"https://plotly.com/r/","desc":"交互式图表库的 R 接口，可生成网页端可交互图形"},{"name":"patchwork","url":"https://patchwork.data-imaginist.com/","desc":"多图拼接工具，用 + 和 / 符号轻松组合 ggplot 图形"},{"name":"RColorBrewer","url":"https://cran.r-project.org/web/packages/RColorBrewer/","desc":"出版级配色方案，提供精心设计的离散色板"}]},{"title":"Python","items":[{"name":"Matplotlib","url":"https://matplotlib.org/","desc":"Python 中广泛使用的绘图库，支持多种图表类型"},{"name":"Seaborn","url":"https://seaborn.pydata.org/","desc":"基于 Matplotlib 的统计可视化库，内置美观主题与统计图"},{"name":"Plotly","url":"https://plotly.com/python/","desc":"交互式可视化库，支持散点、3D、地图等多种图表"},{"name":"Altair","url":"https://altair-viz.github.io/","desc":"基于 Vega-Lite 的声明式统计可视化库"}]}]},{"title":"统计分析","children":[{"title":"R 语言","items":[{"name":"R Stats","url":"https://www.r-project.org/","desc":"R 语言内置的统计分析功能，涵盖广泛的统计方法和模型"},{"name":"rstatix","url":"https://rpkgs.datanovia.com/rstatix/","desc":"管道友好的统计检验封装（t 检验、方差分析、秩和检验等）"},{"name":"broom","url":"https://broom.tidymodels.org/","desc":"把统计模型输出整理为整洁数据框的工具"},{"name":"lme4","url":"https://cran.r-project.org/web/packages/lme4/","desc":"线性与非线性混合效应模型，适用于重复测量数据"}]},{"title":"Python","items":[{"name":"SciPy Stats","url":"https://scipy.org/","desc":"Python 科学计算库 SciPy 中的统计模块，提供多种统计分布和检验方法"},{"name":"statsmodels","url":"https://www.statsmodels.org/","desc":"统计建模库，支持回归、时间序列、假设检验等"},{"name":"pingouin","url":"https://pingouin-stats.org/","desc":"简洁易用的统计检验库，覆盖常用参数与非参数检验"}]}]},{"title":"机器学习","children":[{"title":"Python","items":[{"name":"scikit-learn","url":"https://scikit-learn.org/","desc":"基于 Python 的机器学习库，提供丰富的算法和工具"},{"name":"PyTorch","url":"https://pytorch.org/","desc":"深度学习框架，支持动态计算图和高效的张量运算"},{"name":"TensorFlow","url":"https://www.tensorflow.org/","desc":"深度学习框架，生态成熟，支持生产级部署"},{"name":"XGBoost","url":"https://xgboost.readthedocs.io/","desc":"梯度提升树库，表格数据竞赛与工程的首选"},{"name":"LightGBM","url":"https://lightgbm.readthedocs.io/","desc":"微软出品的梯度提升框架，训练速度快、内存占用低"}]},{"title":"R 语言","items":[{"name":"tidymodels","url":"https://www.tidymodels.org/","desc":"R 语言统一的建模框架，覆盖数据预处理、建模与评估"}]}]}]},{"title":"组学","children":[{"title":"数据平台","children":[{"title":"序列数据库","items":[{"name":"NCBI","url":"https://www.ncbi.nlm.nih.gov/","desc":"美国国家生物技术信息中心，提供 GenBank 等重要数据库"},{"name":"EMBL-EBI","url":"https://www.ebi.ac.uk/","desc":"欧洲分子生物学实验室，提供多种生物信息学资源和工具"},{"name":"DDBJ","url":"https://www.ddbj.nig.ac.jp/","desc":"日本 DNA 数据库，提供核酸、蛋白序列数据的存储和访问"},{"name":"UniProt","url":"https://www.uniprot.org/","desc":"最全面的蛋白质序列与功能注释数据库"},{"name":"Ensembl","url":"https://www.ensembl.org/","desc":"脊椎动物基因组注释与比较基因组学数据库"},{"name":"UCSC Genome Browser","url":"https://genome.ucsc.edu/","desc":"基因组浏览器，支持多物种基因组可视化与自定义轨道"}]},{"title":"测序数据","items":[{"name":"SRA","url":"https://www.ncbi.nlm.nih.gov/sra","desc":"NCBI 原始测序数据仓库，存储高通量测序原始数据"},{"name":"GEO","url":"https://www.ncbi.nlm.nih.gov/geo/","desc":"基因表达数据库，收录芯片与测序表达数据"},{"name":"ENA","url":"https://www.ebi.ac.uk/ena","desc":"欧洲核酸档案库，欧洲的原始测序数据存储中心"}]},{"title":"蛋白互作与通路","items":[{"name":"STRING","url":"https://string-db.org/","desc":"蛋白-蛋白相互作用数据库，支持多种物种"},{"name":"KEGG","url":"https://www.kegg.jp/","desc":"通路、基因与化合物综合数据库"},{"name":"GO","url":"http://geneontology.org/","desc":"基因本体论，提供标准化的基因功能注释体系"}]}]},{"title":"分析软件","children":[{"title":"命令行","items":[{"name":"FastQC","url":"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/","desc":"测序数据质量评估工具，生成可视化质控报告"},{"name":"Trimmomatic","url":"http://www.usadellab.org/cms/?page=trimmomatic","desc":"测序读段修剪工具，去除接头与低质量碱基"},{"name":"cutadapt","url":"https://cutadapt.readthedocs.io/","desc":"命令行接头去除工具，支持多种测序平台"},{"name":"BWA","url":"https://github.com/lh3/bwa","desc":"短读段比对工具，基因组比对的标准选择"},{"name":"STAR","url":"https://github.com/alexdobin/STAR","desc":"RNA-seq 剪接感知比对器，速度极快"},{"name":"HISAT2","url":"https://daehwankimlab.github.io/hisat2/","desc":"基于图索引的 RNA-seq 比对工具"},{"name":"SAMtools","url":"https://www.htslib.org/","desc":"SAM/BAM 文件处理工具集，比对后处理的标准"},{"name":"bcftools","url":"https://samtools.github.io/bcftools/","desc":"VCF/BCF 变异文件处理工具，支持过滤与统计"},{"name":"GATK","url":"https://gatk.broadinstitute.org/","desc":"基因组分析工具包，变异检测的行业标准"},{"name":"MultiQC","url":"https://multiqc.info/","desc":"把多个质控报告汇总为一个交互式 HTML 报告"}]},{"title":"Python","items":[{"name":"Biopython","url":"https://biopython.org/","desc":"生物信息学核心 Python 库，涵盖序列、结构与文件解析"},{"name":"scikit-bio","url":"https://scikit-bio.org/","desc":"用于生物信息学的 Python 库，提供多种生物数据处理和分析工具"},{"name":"scanpy","url":"https://scanpy.readthedocs.io/","desc":"单细胞转录组分析 Python 库，与 AnnData 深度集成"}]},{"title":"R 语言","items":[{"name":"Bioconductor","url":"https://bioconductor.org/","desc":"基于 R 语言的生物信息学软件包集合，涵盖基因组学、转录组学等领域"},{"name":"Seurat","url":"https://satijalab.org/seurat/","desc":"单细胞转录组分析主流 R 包，提供完整分析流程"},{"name":"DESeq2","url":"https://bioconductor.org/packages/release/bioc/html/DESeq2.html","desc":"基于负二项模型的差异表达分析标准工具"},{"name":"edgeR","url":"https://bioconductor.org/packages/release/bioc/html/edgeR.html","desc":"基于经验贝叶斯的差异表达分析工具"},{"name":"limma","url":"https://bioconductor.org/packages/release/bioc/html/limma.html","desc":"线性模型差异分析工具，芯片与测序数据通用"}]}]}]},{"title":"电镜结构解析","children":[{"title":"软件工具","children":[{"title":"三维重构","items":[{"name":"RELION","url":"https://www3.mrc-lmb.cam.ac.uk/groups/scheres/index.php","desc":"基于贝叶斯方法的冷冻电镜三维重构软件，广泛应用于高分辨率结构解析"},{"name":"cryoSPARC","url":"https://cryoem.slac.stanford.edu/","desc":"高性能的冷冻电镜数据处理平台，提供直观的用户界面"},{"name":"cisTEM","url":"https://cistem.org/","desc":"免费开源的冷冻电镜处理软件，支持完整的单颗粒流程"},{"name":"SPHIRE","url":"https://sphire.mpg.de/","desc":"模块化冷冻电镜处理套件，支持高分辨率重构"}]},{"title":"预处理","items":[{"name":"MotionCor2","url":"https://emcore.ucsf.edu/ucsf-software","desc":"逐帧运动校正工具，消除电子束导致的样品漂移"},{"name":"CTFFIND4","url":"https://grigoriefflab.umassmed.edu/ctffind4","desc":"CTF 估计经典工具，评估欠焦与像散"},{"name":"Gctf","url":"https://www.mrc-lmb.cam.ac.uk/kzhang/Gctf/","desc":"GPU 加速的 CTF 估计工具，速度快精度高"},{"name":"Warp","url":"http://www.warpem.com/","desc":"集运动校正、CTF 与颗粒挑选于一体的快速预处理工具"}]},{"title":"颗粒挑选","items":[{"name":"crYOLO","url":"https://cryolo.readthedocs.io/","desc":"基于深度学习的颗粒挑选工具，支持实时挑选"},{"name":"Topaz","url":"https://topaz.readthedocs.io/","desc":"基于 CNN 的颗粒挑选工具，对低信噪比样品表现优异"}]},{"title":"建模与可视化","items":[{"name":"PyMOL","url":"https://pymol.org/","desc":"分子可视化软件，脚本化渲染与结构分析的标准工具"},{"name":"UCSF ChimeraX","url":"https://www.cgl.ucsf.edu/chimerax/","desc":"现代分子可视化软件，密度图处理能力突出"},{"name":"Coot","url":"https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/","desc":"手动模型修正工具，密度图与模型交互调整"},{"name":"phenix","url":"https://phenix-online.org/","desc":"结构精修与验证套件，支持冷冻电镜与晶体学"},{"name":"ISOLDE","url":"https://isolde.cimr.cam.ac.uk/","desc":"交互式实时分子动力学建模工具，基于 ChimeraX"},{"name":"ModelAngelo","url":"https://github.com/3dem/model-angelo","desc":"基于深度学习的密度图自动原子建模工具"}]}]},{"title":"数据库资源","children":[{"title":"结构与密度图","items":[{"name":"EMDB","url":"https://www.ebi.ac.uk/emdb/","desc":"电子显微镜数据银行，存储和分发冷冻电镜三维密度图"},{"name":"EMPIAR","url":"https://www.ebi.ac.uk/empiar/","desc":"冷冻电镜原始数据档案，提供原始电影与颗粒数据"},{"name":"RCSB PDB","url":"https://www.rcsb.org/","desc":"蛋白质数据银行，包含通过X射线、NMR和冷冻电镜确定的生物大分子结构"},{"name":"wwPDB","url":"https://www.wwpdb.org/","desc":"全球 PDB 协调机构，统一管理结构数据标准"},{"name":"PDBe","url":"https://www.ebi.ac.uk/pdbe/","desc":"欧洲 PDB 节点，提供结构数据与可视化工具"},{"name":"AlphaFold DB","url":"https://alphafold.ebi.ac.uk/","desc":"AlphaFold 预测蛋白质结构数据库，覆盖数亿蛋白质"}]}]}]}]'),Ob=Object.freeze(Object.defineProperty({__proto__:null,default:Db},Symbol.toStringTag,{value:"Module"})),Ib=[{name:"Advanced",count:1},{name:"AlphaFold",count:1},{name:"Atomic Modeling",count:1},{name:"AutoDock Vina",count:1},{name:"Bash",count:1},{name:"ChimeraX",count:1},{name:"Command Line",count:1},{name:"Computational Biology",count:2},{name:"Coot",count:1},{name:"Cryo-Electron Microscopy",count:2},{name:"cryo-EM",count:2},{name:"Data Processing",count:2},{name:"Data Visualization",count:1},{name:"dplyr",count:1},{name:"Getting Started",count:2},{name:"ggplot2",count:1},{name:"Linux",count:1},{name:"Molecular Docking",count:1},{name:"NumPy",count:1},{name:"Pandas",count:1},{name:"phenix",count:1},{name:"Protein Design",count:1},{name:"PyMOL",count:1},{name:"Python",count:3},{name:"R Language",count:3},{name:"RELION",count:1},{name:"Review",count:1},{name:"RFdiffusion",count:1},{name:"Rosetta",count:1},{name:"Scripting",count:1},{name:"Shell",count:1},{name:"Structure Refinement",count:1},{name:"Structure Visualization",count:1},{name:"tidyverse",count:1},{name:"Tutorial",count:10},{name:"Virtual Screening",count:1}],Nb=Object.freeze(Object.defineProperty({__proto__:null,default:Ib},Symbol.toStringTag,{value:"Module"})),Fb=[{name:"蛋白质设计",count:1},{name:"分子对接",count:1},{name:"计算生物学",count:2},{name:"脚本编程",count:1},{name:"教程",count:10},{name:"结构精修",count:1},{name:"结构可视化",count:1},{name:"进阶",count:1},{name:"冷冻电镜",count:2},{name:"命令行",count:1},{name:"入门",count:2},{name:"数据处理",count:2},{name:"数据可视化",count:1},{name:"虚拟筛选",count:1},{name:"原子建模",count:1},{name:"综述",count:1},{name:"AlphaFold",count:1},{name:"AutoDock Vina",count:1},{name:"Bash",count:1},{name:"ChimeraX",count:1},{name:"Coot",count:1},{name:"cryo-EM",count:2},{name:"dplyr",count:1},{name:"ggplot2",count:1},{name:"Linux",count:1},{name:"NumPy",count:1},{name:"Pandas",count:1},{name:"phenix",count:1},{name:"PyMOL",count:1},{name:"Python",count:3},{name:"R语言",count:3},{name:"RELION",count:1},{name:"RFdiffusion",count:1},{name:"Rosetta",count:1},{name:"Shell",count:1},{name:"tidyverse",count:1}],$b=Object.freeze(Object.defineProperty({__proto__:null,default:Fb},Symbol.toStringTag,{value:"Module"})),Bb=[{name:"StructureOfProteinDemo",title:"Protein Structure Determination Demo Project",desc:"Protein structure placeholder demo",date:"2026-04-12",tags:["Protein","Structure","Cryo-EM"],doi:"10.1234/demo-structure.2026",url:"https://www.demo-structure.org",status:"completed",journal:"Journal of Demo Science",year:2026,authors:["Demo Author A","Demo Author B"]}],qb=Object.freeze(Object.defineProperty({__proto__:null,default:Bb},Symbol.toStringTag,{value:"Module"})),zb=[{name:"StructureOfProteinDemo",title:"蛋白质结构解析演示课题",desc:"蛋白质结构占位demo",date:"2026-04-12",tags:["Protein","Structure","Cryo-EM"],doi:"10.1234/demo-structure.2026",url:"https://www.demo-structure.org",status:"completed",journal:"Journal of Demo Science",year:2026,authors:["Demo Author A","Demo Author B"]}],Vb=Object.freeze(Object.defineProperty({__proto__:null,default:zb},Symbol.toStringTag,{value:"Module"})),Ub=`<h1>Mainstream Tools for Protein Design</h1>
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
<p>The next article will introduce mainstream tools for virtual screening: a complete toolchain for docking small molecules to targets.</p>`,Hb=`<h1>蛋白质设计主流工具</h1>
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
<p>下一篇将介绍虚拟筛选的主流工具：把小分子与靶点对接的完整工具链。</p>`,Gb=`<h1>Mainstream Tools for Virtual Screening</h1>
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
<p>This completes the computational biology direction: protein design tools + virtual screening tools, with two main threads forming a closed loop of "design-screen" computational drug discovery.</p>`,Wb=`<h1>虚拟筛选主流工具</h1>
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
<p>至此计算生物学方向完成：蛋白质设计工具 + 虚拟筛选工具，两条主线构成"设计-筛选"的计算药物发现闭环。</p>`,Kb=`<h1>Bash Programming: Variables, Conditionals, Loops, and Practical Scripts</h1>
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
<p>At this point, the programming track tutorials are complete: Python ×3, R ×3, Linux ×1, Bash ×1, from zero to hands-on practice.</p>`,Xb=`<h1>Bash 编程：变量、条件、循环与实用脚本</h1>
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
<p>至此编程方向教程完成：Python ×3、R ×3、Linux ×1、Bash ×1，从零到实战。</p>`,Yb=`<h1>Linux Command Line Basics: File System, Permissions, and Text Processing</h1>
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
<p>The next article will introduce Bash programming: organizing commands into reusable scripts.</p>`,Qb=`<h1>Linux 命令行基础：文件系统、权限与文本处理</h1>
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
<p>下一篇将介绍 Bash 编程：把命令组织成可复用的脚本。</p>`,Jb=`<h1>Python Advanced: Functions, Classes, and Modules</h1>
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
<p>The next article will cover Python data processing in practice: file I/O, regular expressions, and NumPy/Pandas.</p>`,Zb=`<h1>Python 进阶：函数、类与模块</h1>
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
<p>下一篇将介绍 Python 数据处理实战：文件 IO、正则表达式与 NumPy/Pandas。</p>`,s_=`<h1>Introduction to Python Programming: Environment, Syntax, and Data Types</h1>
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
<p>The next article will introduce functions, classes, and modules, moving into real engineering-style programming.</p>`,n_=`<h1>Python 编程入门：环境、语法与数据类型</h1>
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
<p>下一篇将介绍函数、类与模块，进入真正的工程化编程。</p>`,a_=`<h1>Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas</h1>
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
<p>At this point, the Python tutorial trilogy is complete: Beginner → Intermediate → Data Processing.</p>`,e_=`<h1>Python 数据处理实战：文件 IO、正则与 NumPy/Pandas</h1>
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
<p>至此 Python 教程三部曲完成：入门 → 进阶 → 数据处理。</p>`,t_=`<h1>R Language Primer: Data Structures, Vectorization, and Functions</h1>
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
<p>The next article will introduce tidyverse: using <code>dplyr</code> for elegant data manipulation.</p>`,l_=`<h1>R 语言入门：数据结构、向量化与函数</h1>
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
<p>下一篇将介绍 tidyverse：用 <code>dplyr</code> 优雅地完成数据操作。</p>`,o_=`<h1>R ggplot2 Data Visualization: Grammar of Graphics and Common Charts</h1>
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
<p>With this, the R tutorial trilogy is complete: Introduction → tidyverse → ggplot2.</p>`,p_=`<h1>R ggplot2 数据可视化：图层语法与常用图表</h1>
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
<p>至此 R 教程三部曲完成：入门 → tidyverse → ggplot2。</p>`,c_=`<h1>R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning</h1>
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
<p>The next article will introduce ggplot2: creating publication-quality graphics using the grammar of layers.</p>`,r_=`<h1>R tidyverse 数据操作：dplyr 管道与数据清洗</h1>
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
<p>下一篇将介绍 ggplot2：用图层语法绘制出版级图表。</p>`,i_=`<h1>Review of Cryo-EM Technology</h1>
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
<p>The past decade of cryo-EM is a model of synergy among engineering, physics, and computational biology. In the next decade, it will join forces with AI to redefine how we understand life.</p>`,u_=`<h1>冷冻电镜技术综述</h1>
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
<p>冷冻电镜的十年，是工程、物理与计算生物学协同的典范。下一个十年，它将与 AI 一起重新定义我们理解生命的方式。</p>`,h_=`<h1>Cryo-EM Single Particle Analysis: Full Data Processing Workflow</h1>
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
<p>Future articles will cover structure visualization (PyMOL/ChimeraX) and atomic modeling (Coot/phenix), turning density maps into atomic models.</p>`,d_=`<h1>冷冻电镜单颗粒分析：数据处理全流程</h1>
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
<p>后续文章将介绍结构可视化（PyMOL/ChimeraX）与原子建模（Coot/phenix），把密度图变成原子模型。</p>`,m_=`<h1>Atomic Modeling and Refinement</h1>
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
<p>This completes the structural biology pipeline: data processing → technical review → visualization → atomic modeling, forming a complete loop from experiment to model.</p>`,j_=`<h1>原子建模与精修</h1>
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
<p>至此结构生物学方向完成：数据处理流程 → 技术综述 → 可视化 → 原子建模，构成从实验到模型的完整闭环。</p>`,g_=`<h1>Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice</h1>
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
<p>The next article will cover atomic modeling: how to build and refine atomic models from density maps.</p>`,f_=`<h1>生物大分子结构可视化：PyMOL 与 ChimeraX 实战</h1>
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
<p>下一篇将介绍原子建模：如何从密度图搭建并精修原子模型。</p>`,Ht=()=>{const s=globalThis.__GBLOG_LOCALE__;return s||(typeof window<"u"?localStorage.getItem("locale"):null)||"zh-CN"},b_=(s,n=".json")=>{const e=Ht()==="en-US";return s.endsWith("-en")?`${s}${n}`:e?`${s}-en${n}`:`${s}${n}`},Ze=Object.assign({"/data-branch/content/about-en.json":db,"/data-branch/content/about.json":jb,"/data-branch/content/categories-en.json":fb,"/data-branch/content/categories.json":_b,"/data-branch/content/notes-en.json":vb,"/data-branch/content/notes.json":Cb,"/data-branch/content/posts-en.json":xb,"/data-branch/content/posts.json":Pb,"/data-branch/content/projects-en.json":Tb,"/data-branch/content/projects.json":Rb,"/data-branch/content/resources-en.json":Mb,"/data-branch/content/resources.json":Ob,"/data-branch/content/tags-en.json":Nb,"/data-branch/content/tags.json":$b,"/data-branch/content/topics-en.json":qb,"/data-branch/content/topics.json":Vb}),Ac=Object.assign({"/data-branch/content/html/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en.html":Ub,"/data-branch/content/html/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools.html":Hb,"/data-branch/content/html/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools-en.html":Gb,"/data-branch/content/html/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools.html":Wb,"/data-branch/content/html/notes/Programming/bash/bash-scripting/bash-scripting-en.html":Kb,"/data-branch/content/html/notes/Programming/bash/bash-scripting/bash-scripting.html":Xb,"/data-branch/content/html/notes/Programming/linux/linux-basics/linux-basics-en.html":Yb,"/data-branch/content/html/notes/Programming/linux/linux-basics/linux-basics.html":Qb,"/data-branch/content/html/notes/Programming/python/python-advanced/python-advanced-en.html":Jb,"/data-branch/content/html/notes/Programming/python/python-advanced/python-advanced.html":Zb,"/data-branch/content/html/notes/Programming/python/python-basics/python-basics-en.html":s_,"/data-branch/content/html/notes/Programming/python/python-basics/python-basics.html":n_,"/data-branch/content/html/notes/Programming/python/python-data/python-data-en.html":a_,"/data-branch/content/html/notes/Programming/python/python-data/python-data.html":e_,"/data-branch/content/html/notes/Programming/r/r-basics/r-basics-en.html":t_,"/data-branch/content/html/notes/Programming/r/r-basics/r-basics.html":l_,"/data-branch/content/html/notes/Programming/r/r-ggplot2/r-ggplot2-en.html":o_,"/data-branch/content/html/notes/Programming/r/r-ggplot2/r-ggplot2.html":p_,"/data-branch/content/html/notes/Programming/r/r-tidyverse/r-tidyverse-en.html":c_,"/data-branch/content/html/notes/Programming/r/r-tidyverse/r-tidyverse.html":r_,"/data-branch/content/html/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en.html":i_,"/data-branch/content/html/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview.html":u_,"/data-branch/content/html/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow-en.html":h_,"/data-branch/content/html/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow.html":d_,"/data-branch/content/html/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling-en.html":m_,"/data-branch/content/html/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling.html":j_,"/data-branch/content/html/notes/StructuralBiology/visualization/structure-visualization/structure-visualization-en.html":g_,"/data-branch/content/html/notes/StructuralBiology/visualization/structure-visualization/structure-visualization.html":f_}),Rc=Object.assign({"/data-branch/cache/en/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en.md":()=>fs(()=>import("./protein-design-tools-en-kwRRVDCZ.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools-en.md":()=>fs(()=>import("./virtual-screening-tools-en-Ds1a_3n3.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/bash/bash-scripting/bash-scripting-en.md":()=>fs(()=>import("./bash-scripting-en-DWq6sXnJ.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/linux/linux-basics/linux-basics-en.md":()=>fs(()=>import("./linux-basics-en-DSb4s4PH.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/python/python-advanced/python-advanced-en.md":()=>fs(()=>import("./python-advanced-en-BQKRq7MH.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/python/python-basics/python-basics-en.md":()=>fs(()=>import("./python-basics-en-BvTX6dVV.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/python/python-data/python-data-en.md":()=>fs(()=>import("./python-data-en-DTE7E1vm.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/r/r-basics/r-basics-en.md":()=>fs(()=>import("./r-basics-en-Dzr6Fzbi.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/r/r-ggplot2/r-ggplot2-en.md":()=>fs(()=>import("./r-ggplot2-en-CvP2vjkk.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/r/r-tidyverse/r-tidyverse-en.md":()=>fs(()=>import("./r-tidyverse-en-CuLVIm9c.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview-en.md":()=>fs(()=>import("./cryoem-overview-en-DsP4yYlt.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow-en.md":()=>fs(()=>import("./cryoem-workflow-en-BpFgioRH.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling-en.md":()=>fs(()=>import("./atomic-modeling-en-BX1y6yZD.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/StructuralBiology/visualization/structure-visualization/structure-visualization-en.md":()=>fs(()=>import("./structure-visualization-en-CHXUkj4Z.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools-en.md":()=>fs(()=>import("./protein-design-tools-en-BdlMbrNY.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/ComputationalBiology/protein-design/protein-design-tools/protein-design-tools.md":()=>fs(()=>import("./protein-design-tools-Jzs48-a0.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools-en.md":()=>fs(()=>import("./virtual-screening-tools-en-BEtVbCe5.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/ComputationalBiology/virtual-screening/virtual-screening-tools/virtual-screening-tools.md":()=>fs(()=>import("./virtual-screening-tools-BdHWOXdu.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/bash/bash-scripting/bash-scripting-en.md":()=>fs(()=>import("./bash-scripting-en-CYovvjQL.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/bash/bash-scripting/bash-scripting.md":()=>fs(()=>import("./bash-scripting-DEkzGDwa.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/linux/linux-basics/linux-basics-en.md":()=>fs(()=>import("./linux-basics-en-BDL0HlS2.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/linux/linux-basics/linux-basics.md":()=>fs(()=>import("./linux-basics-UyN4mB_u.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-advanced/python-advanced-en.md":()=>fs(()=>import("./python-advanced-en-D3xdyM1A.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-advanced/python-advanced.md":()=>fs(()=>import("./python-advanced-iN1UORRl.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-basics/python-basics-en.md":()=>fs(()=>import("./python-basics-en-DOB3GCxE.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-basics/python-basics.md":()=>fs(()=>import("./python-basics-uh9Cz9FQ.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-data/python-data-en.md":()=>fs(()=>import("./python-data-en-BH1w0rP8.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-data/python-data.md":()=>fs(()=>import("./python-data-CmFq1YcA.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/r/r-basics/r-basics-en.md":()=>fs(()=>import("./r-basics-en-CPaztwoM.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/r/r-basics/r-basics.md":()=>fs(()=>import("./r-basics-CqBb8X8p.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/r/r-ggplot2/r-ggplot2.md":()=>fs(()=>import("./r-ggplot2-CK3q5J49.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/r/r-tidyverse/r-tidyverse-en.md":()=>fs(()=>import("./r-tidyverse-en-C13YYJSL.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/r/r-tidyverse/r-tidyverse.md":()=>fs(()=>import("./r-tidyverse-v7eXFCz-.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/StructuralBiology/cryoem/cryoem-overview/cryoem-overview.md":()=>fs(()=>import("./cryoem-overview-DT2Crx5Y.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/StructuralBiology/cryoem/cryoem-workflow/cryoem-workflow.md":()=>fs(()=>import("./cryoem-workflow-C2zs9nW6.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/StructuralBiology/modeling/atomic-modeling/atomic-modeling.md":()=>fs(()=>import("./atomic-modeling-BGRUqNS4.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/StructuralBiology/visualization/structure-visualization/structure-visualization.md":()=>fs(()=>import("./structure-visualization-B3-STJgI.js"),[]).then(s=>s.default)}),se=(s,n)=>{const a=b_(s,".json"),e=Object.keys(Ze).find(t=>t.includes(a));if(e)return Ze[e].default;if(console.error(`Failed to load JSON content: ${a}`),Ht()==="en-US"){const t=Object.keys(Ze).find(l=>l.includes(`${s}.json`));if(t)return Ze[t].default}return n},__=s=>{const a=Ht()==="en-US",e=s.replace(/\.md$/i,""),t=[];e.endsWith("-en")?(t.push(e),a&&t.push(e.replace(/-en$/,""))):(t.push(a?`${e}-en`:e),t.push(e));for(const l of t){const o=Object.keys(Ac).find(p=>p.includes(`/content/html/${l}.html`));if(o)return Ac[o]}return`<h1>Content Not Available</h1>
<p>The requested content could not be loaded.</p>`},Gt=()=>se("categories",[]),y_=()=>se("posts",[]),v_=()=>se("notes",[]),w_=()=>se("tags",[]),C_=()=>se("resources",[]),iu={introduction:"",experience:[],section:[],contacts:[]},k_=()=>se("about",iu),x_=async s=>{const n=Ht()==="en-US"?"-en":"",a=`${s.replace(/\.md$/,"")}${n}.md`,e=Object.keys(Rc).find(t=>t.endsWith(`/${a}`)||t.endsWith(a));if(!e)return null;try{return await Rc[e]()}catch(t){return console.error(`Failed to load markdown source: ${a}`,t),null}},S_={class:"navigation-tree"},P_={class:"tree-label"},E_={class:"article-list"},T_=["aria-expanded","onClick"],A_={class:"tree-item__text"},R_={class:"article-list tree-sublist"},L_=["aria-expanded","onClick"],M_={class:"tree-item__text"},D_={class:"article-list tree-sublist"},O_=Ns({__name:"NavigationTree",setup(s){const n=ya(),a=rs([]),e=rs(""),t=rs([]),{data:l}=fa(()=>Gt(),[]);function o(h){return Jn(nu(h))}function p(h){const d=e.value.replace(/\.md$/,"").replace(/-en$/,""),b=h.replace(/\.md$/,"").replace(/-en$/,"");return d===b}function c(h){return!t.value.includes(h.id)}function u(h){t.value=c(h)?[...t.value,h.id]:t.value.filter(d=>d!==h.id)}function r(){const h=e.value;if(!h){a.value=[];return}const d=h.split("/").filter(Boolean);if(d.length<2){a.value=[];return}const b=d[0],f=d[1];let C=null;s:for(const x of l.value)for(const R of x.items){if(R.name!==f)continue;const E=P=>P.articleUrl.includes(`/article/${b}/`);if(R.articles?.some(E)||R.categories.some(P=>P.articles.some(U=>U.articleUrl.includes(`/article/${b}/${f}/`)))){C=R;break s}}if(!C){a.value=[];return}const m=[],g=[],j=(x,R)=>({title:x,path:Xl(R)});C.articles?.forEach(x=>{x.articleUrl&&m.push(j(x.title,x.articleUrl))}),C.categories.forEach(x=>{const R=[];x.articles.forEach(E=>{E.articleUrl&&R.push(j(E.title,E.articleUrl))}),R.length&&g.push({id:`${C.name}/${x.key||x.title}`,name:x.title||x.key,files:R})});const y=[{name:C.title||C.name||f,files:m,children:g}];a.value=y}function i(){e.value=ct(n.params.path),r()}return Vs(()=>n.params.path,i,{immediate:!0}),Vs(l,()=>r()),En(()=>{r()}),(h,d)=>{const b=Za("router-link");return I(),N("div",S_,[(I(!0),N(gs,null,As(a.value,f=>(I(),N("div",{key:f.name,class:"tree-group"},[S("div",P_,G(f.name),1),S("ul",E_,[(I(!0),N(gs,null,As(f.children,C=>(I(),N("li",{key:C.id,class:"article-item"},[S("button",{type:"button",class:Ts(["tree-item tree-item--folder",{"tree-item--branch-open":c(C)}]),"aria-expanded":c(C),onClick:m=>u(C)},[S("span",A_,G(C.name),1),S("i",{class:Ts(["fas fa-chevron-right tree-caret",{"tree-caret--open":c(C)}]),"aria-hidden":"true"},null,2)],10,T_),an(S("ul",R_,[(I(!0),N(gs,null,As(C.files,m=>(I(),N("li",{key:m.path,class:"article-item"},[bs(b,{to:o(m.path),class:Ts(["tree-item tree-item--l2",{"tree-item--active":p(m.path)}])},{default:jn(()=>[Cs(G(m.title),1)]),_:2},1032,["to","class"])]))),128)),(I(!0),N(gs,null,As(C.children,m=>(I(),N("li",{key:m.id,class:"article-item"},[S("button",{type:"button",class:Ts(["tree-item tree-item--folder tree-item--l2",{"tree-item--branch-open":c(m)}]),"aria-expanded":c(m),onClick:g=>u(m)},[S("span",M_,G(m.name),1),S("i",{class:Ts(["fas fa-chevron-right tree-caret",{"tree-caret--open":c(m)}]),"aria-hidden":"true"},null,2)],10,L_),an(S("ul",D_,[(I(!0),N(gs,null,As(m.files,g=>(I(),N("li",{key:g.path,class:"article-item"},[bs(b,{to:o(g.path),class:Ts(["tree-item tree-item--l3",{"tree-item--active":p(g.path)}])},{default:jn(()=>[Cs(G(g.title),1)]),_:2},1032,["to","class"])]))),128))],512),[[yt,c(m)]])]))),128))],512),[[yt,c(C)]])]))),128)),(I(!0),N(gs,null,As(f.files,C=>(I(),N("li",{key:C.path,class:"article-item"},[bs(b,{to:o(C.path),class:Ts(["tree-item tree-item--l1",{"tree-item--active":p(C.path)}])},{default:jn(()=>[Cs(G(C.title),1)]),_:2},1032,["to","class"])]))),128))])]))),128))])}}}),hn=(s,n)=>{const a=s.__vccOpts||s;for(const[e,t]of n)a[e]=t;return a},uu=hn(O_,[["__scopeId","data-v-c348384c"]]),I_=["aria-label"],N_=["aria-label"],F_={key:0,class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},$_={key:1,class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},B_=["aria-label"],q_=["aria-label"],z_=Ns({__name:"NavActions",props:{mobile:{type:Boolean,default:!1}},emits:["open-search","toggle-menu"],setup(s,{emit:n}){const a=n,{t:e}=Us(),t=ya(),l=qe(),o=Do(),p=ns(()=>o.theme==="dark"?!0:o.theme==="light"?!1:typeof window<"u"&&window.matchMedia("(prefers-color-scheme: dark)").matches),c=i=>{i.currentTarget?.blur()},u=()=>o.toggleTheme(),r=()=>ub(l,t);return(i,h)=>(I(),N("div",{class:Ts(["d-flex",s.mobile?"d-lg-none ms-auto app-nav__actions app-nav__actions--mobile":"d-none d-lg-flex ms-auto app-nav__actions"])},[S("button",{class:"icon-btn",onClick:h[0]||(h[0]=d=>a("open-search")),onFocus:c,"aria-label":ss(e)("search")},[...h[2]||(h[2]=[S("svg",{class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[S("circle",{cx:"11",cy:"11",r:"8"}),S("line",{x1:"21",y1:"21",x2:"16.65",y2:"16.65"})],-1)])],40,I_),S("button",{class:"icon-btn",onClick:u,onFocus:c,"aria-label":ss(e)("theme")},[p.value?(I(),N("svg",F_,[...h[3]||(h[3]=[tp('<circle cx="12" cy="12" r="5" data-v-da55aa12></circle><line x1="12" y1="1" x2="12" y2="3" data-v-da55aa12></line><line x1="12" y1="21" x2="12" y2="23" data-v-da55aa12></line><line x1="4.22" y1="4.22" x2="5.64" y2="5.64" data-v-da55aa12></line><line x1="18.36" y1="18.36" x2="19.78" y2="19.78" data-v-da55aa12></line><line x1="1" y1="12" x2="3" y2="12" data-v-da55aa12></line><line x1="21" y1="12" x2="23" y2="12" data-v-da55aa12></line><line x1="4.22" y1="19.78" x2="5.64" y2="18.36" data-v-da55aa12></line><line x1="18.36" y1="5.64" x2="19.78" y2="4.22" data-v-da55aa12></line>',9)])])):(I(),N("svg",$_,[...h[4]||(h[4]=[S("path",{d:"M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"},null,-1)])]))],40,N_),S("button",{class:"icon-btn",onClick:r,onFocus:c,"aria-label":ss(e)("language")},[...h[5]||(h[5]=[tp('<svg class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.75" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" data-v-da55aa12><path d="m5 8 6 6" data-v-da55aa12></path><path d="m4 14 6-6 2-3" data-v-da55aa12></path><path d="M2 5h12" data-v-da55aa12></path><path d="M7 2h1" data-v-da55aa12></path><path d="m22 22-5-10-5 10" data-v-da55aa12></path><path d="M14 18h6" data-v-da55aa12></path></svg>',1)])],40,B_),s.mobile?(I(),N("button",{key:0,class:"icon-btn",onClick:h[1]||(h[1]=d=>a("toggle-menu")),onFocus:c,"aria-label":ss(e)("menu")},[...h[6]||(h[6]=[S("svg",{class:"app-nav__icon",width:"18",height:"18",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[S("line",{x1:"3",y1:"12",x2:"21",y2:"12"}),S("line",{x1:"3",y1:"6",x2:"21",y2:"6"}),S("line",{x1:"3",y1:"18",x2:"21",y2:"18"})],-1)])],40,q_)):ms("",!0)],2))}}),Lc=hn(z_,[["__scopeId","data-v-da55aa12"]]),V_={class:"app-header"},U_={class:"container px-0"},H_={class:"navbar navbar-expand-lg app-nav"},G_={class:"container-fluid d-flex app-nav__inner"},W_={class:"app-nav__wordmark"},K_={class:"navbar-nav mb-2 mb-lg-0 app-nav__links"},X_={key:0,class:"app-nav__link-divider","aria-hidden":"true"},Y_={class:"offcanvas-panel"},Q_={class:"offcanvas-section"},J_={class:"offcanvas-head"},Z_={class:"offcanvas-card"},sy={class:"list-unstyled m-0"},ny={key:0,class:"offcanvas-section"},ay={class:"offcanvas-head"},ey=Ns({__name:"AppHeader",emits:["open-search"],setup(s,{emit:n}){const a=ya(),{t:e}=Us(),t=Do(),l=Mo(),o=rs(!1),p=rs(!1),c=n,u=ns(()=>[{text:e("categories"),href:Jn("/category")},{text:e("resources"),href:Jn("/resource")},{text:e("about"),href:Jn("/about")}]),r=C=>C===a.path?!0:C!==Jn("/")&&a.path.startsWith(C),i=ns(()=>a.path.includes("/article/")),h=()=>{window.innerWidth<992?d():o.value=!o.value},d=()=>{rb(),p.value=!0},b=()=>{p.value=!1,Ec()},f=C=>{const m=C.target;m&&m.closest("a")&&b()};return En(()=>{t.initTheme(),l.initLocale()}),zn(()=>{Ec()}),(C,m)=>{const g=Za("RouterLink");return I(),N(gs,null,[S("header",V_,[S("div",U_,[S("nav",H_,[S("div",G_,[bs(g,{class:"navbar-brand app-nav__brand",to:ss(Jn)("/"),onClick:m[0]||(m[0]=j=>o.value=!1)},{default:jn(()=>[S("span",W_,[Cs(G(ss(cn).author),1),m[4]||(m[4]=S("span",{class:"app-nav__apos"},"’",-1)),m[5]||(m[5]=Cs("s blog",-1))])]),_:1},8,["to"]),bs(Lc,{mobile:"",onOpenSearch:m[1]||(m[1]=j=>c("open-search")),onToggleMenu:h}),S("div",{class:Ts(["navbar-collapse collapse",{show:o.value}])},[S("ul",K_,[u.value.length?(I(),N("li",X_)):ms("",!0),(I(!0),N(gs,null,As(u.value,j=>(I(),N("li",{class:"nav-item",key:j.text},[bs(g,{class:Ts(["nav-link app-nav__link",{active:r(j.href)}]),to:j.href,onClick:m[2]||(m[2]=y=>o.value=!1)},{default:jn(()=>[Cs(G(j.text),1)]),_:2},1032,["to","class"])]))),128))])],2),bs(Lc,{onOpenSearch:m[3]||(m[3]=j=>c("open-search"))})])])])]),p.value?(I(),N("div",{key:0,class:"mobile-offcanvas d-lg-none",onClick:vn(b,["self"])},[S("div",Y_,[S("div",Q_,[S("div",J_,G(ss(e)("menu")),1),S("div",Z_,[S("ul",sy,[(I(!0),N(gs,null,As(u.value,j=>(I(),N("li",{key:j.text,class:"my-1"},[bs(g,{to:j.href,class:Ts(["offcanvas-link d-flex align-items-center",{active:r(j.href)}]),onClick:b},{default:jn(()=>[S("span",null,G(j.text),1),m[6]||(m[6]=S("i",{class:"fas fa-chevron-right offcanvas-link__chevron"},null,-1))]),_:2},1032,["to","class"])]))),128))])])]),i.value?(I(),N("div",ny,[S("div",ay,G(ss(e)("tableOfContents")),1),S("div",{class:"offcanvas-tree offcanvas-card",onClick:f},[bs(uu)])])):ms("",!0)]),S("div",{class:"offcanvas-backdrop",onClick:b})])):ms("",!0)],64)}}}),ty=hn(ey,[["__scopeId","data-v-b943a66d"]]),ly={class:"site-footer"},oy={class:"container"},py={class:"site-footer__inner"},cy={class:"footer-copy"},ry={class:"footer-designed"},iy={class:"footer-designed__text"},uy={class:"footer-designed__name"},hy={class:"footer-designed__icons"},dy=["href"],my=["href"],jy=Ns({__name:"AppFooter",setup(s){const{t:n}=Us(),a=new Date().getFullYear(),e=cn.startYear&&cn.startYear<a?`${cn.startYear} - ${a}`:`${a}`;return(t,l)=>(I(),N("footer",ly,[S("div",oy,[S("div",py,[S("p",cy,"© "+G(ss(e))+" "+G(ss(cn).author),1),S("div",ry,[S("p",iy,[Cs(G(ss(n)("designedByPrefix")),1),S("strong",uy,G(ss(cn).author),1),Cs(G(ss(n)("designedBySuffix")),1)]),S("div",hy,[S("a",{href:ss(cn).github,target:"_blank",rel:"noopener noreferrer",class:"footer-link","aria-label":"GitHub"},[...l[0]||(l[0]=[S("svg",{class:"footer-icon",width:"16",height:"16",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[S("path",{d:"M9 19c-5 1.5-5-2.5-7-3m14 6v-3.87a3.37 3.37 0 0 0-.94-2.61c3.14-.35 6.44-1.54 6.44-7A5.44 5.44 0 0 0 20 4.77 5.07 5.07 0 0 0 19.91 1S18.73.65 16 2.48a13.38 13.38 0 0 0-7 0C6.27.65 5.09 1 5.09 1A5.07 5.07 0 0 0 5 4.77a5.44 5.44 0 0 0-1.5 3.78c0 5.42 3.3 6.61 6.44 7A3.37 3.37 0 0 0 9 18.13V22"})],-1)])],8,dy),S("a",{href:`mailto:${ss(cn).email}`,class:"footer-link","aria-label":"Email"},[...l[1]||(l[1]=[S("svg",{class:"footer-icon",width:"16",height:"16",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true"},[S("path",{d:"M4 4h16c1.1 0 2 .9 2 2v12c0 1.1-.9 2-2 2H4c-1.1 0-2-.9-2-2V6c0-1.1.9-2 2-2z"}),S("polyline",{points:"22,6 12,13 2,6"})],-1)])],8,my)])])])])]))}}),gy=hn(jy,[["__scopeId","data-v-443df76b"]]);function hu(s){const{sourceId:n,mode:a,onRelease:e}=s,t=s.buttonHeight??40,l=s.gap??48,o=s.margin??20,p=rs(!1),c=rs(null),u=rs(!1),r=rs(0),i=rs(0),h=rs(s.defaultTop),d=rs(!1),b=()=>({gap:l,minTop:a==="stack"?o:o+l,maxTop:a==="stack"?Math.max(0,window.innerHeight-t-o-l):window.innerHeight-t-o}),f=P=>{const{minTop:U,maxTop:Y}=b();return Math.max(U,Math.min(Y,P))},C=P=>a==="stack"?P+l:P,m=P=>a==="stack"?P-l:P;function g(){c.value=C(h.value),!p.value&&(p.value=!0,requestAnimationFrame(()=>{c.value!==null&&window.dispatchEvent(new CustomEvent("floating-buttons-base-top",{detail:{baseTop:c.value,source:n}})),p.value=!1}))}function j(P){const U=P.detail,Y=U?.baseTop;U?.source!==n&&typeof Y=="number"&&(h.value=f(m(Y)))}function y(P){P.preventDefault(),u.value=!0,d.value=!1,r.value=P.touches[0].clientY,i.value=h.value}function x(P){d.value=!0;const U=P.touches[0].clientY-r.value;h.value=f(i.value+U),g(),P.preventDefault()}function R(P){P.preventDefault(),u.value=!1,d.value||e()}function E(){window.addEventListener("floating-buttons-base-top",j)}function O(){window.removeEventListener("floating-buttons-base-top",j)}return{isDragging:u,buttonTop:h,clampTop:f,dispatchBaseTop:g,onTouchStart:y,onTouchMove:x,onTouchEnd:R,subscribe:E,unsubscribe:O}}const fy=["aria-label"],by=Ns({__name:"BackToTop",setup(s){const{t:n}=Us(),a=rs(!1),e=rs(null),t=rs(null);function l(){const j=typeof window.matchMedia=="function"&&window.matchMedia("(prefers-reduced-motion: reduce)").matches;window.scrollTo({top:0,behavior:j?"auto":"smooth"})}const{isDragging:o,buttonTop:p,dispatchBaseTop:c,onTouchStart:u,onTouchMove:r,onTouchEnd:i,subscribe:h,unsubscribe:d}=hu({sourceId:"btt",defaultTop:(typeof window<"u"?window.innerHeight:1024)-100,mode:"match",onRelease:l}),b=()=>{o.value||l()},f=u,C=r,m=i;function g(){if(e.value)return;const j=document.createElement("div");j.style.cssText="position:absolute;left:0;top:0;width:1px;height:181px;pointer-events:none;visibility:hidden",document.body.appendChild(j),t.value=j,e.value=new IntersectionObserver(([y])=>{a.value=!y.isIntersecting,a.value&&c()},{threshold:0}),e.value.observe(j)}return En(()=>{g(),h(),p.value=window.innerHeight-100,c()}),zn(()=>{e.value&&(e.value.disconnect(),e.value=null),t.value&&(t.value.remove(),t.value=null),d()}),(j,y)=>an((I(),N("button",{class:"back-to-top d-flex align-items-center justify-content-center",onClick:b,"aria-label":ss(n)("backToTop"),onTouchstart:y[0]||(y[0]=vn((...x)=>ss(f)&&ss(f)(...x),["prevent","stop"])),onTouchmove:y[1]||(y[1]=vn((...x)=>ss(C)&&ss(C)(...x),["prevent","stop"])),onTouchend:y[2]||(y[2]=vn((...x)=>ss(m)&&ss(m)(...x),["prevent","stop"])),style:qn({top:ss(p)+"px"})},[...y[3]||(y[3]=[S("i",{class:"fas fa-arrow-up"},null,-1)])],44,fy)),[[yt,a.value]])}}),_y=hn(by,[["__scopeId","data-v-0f89b0e5"]]),yy={class:"page-wrap"},vy={class:"main-content"},wy=Ns({__name:"App",setup(s){const n=rs(!1),a=dh(()=>fs(()=>import("./SearchModal-BsN3_RSp.js"),__vite__mapDeps([0,1]))),e=ya(),t=Mo(),{locale:l}=Us(),o=ns(()=>{const i=Zi(e.path);return i==="/"||i===""?"":i}),p=ns(()=>Vt(ea(l.value)));It({htmlAttrs:{lang:()=>xt[ea(l.value)]},link:ns(()=>{const i=o.value.includes("/article/"),h=i?o.value.replace(/-en$/,""):o.value,d=i&&!o.value.endsWith("-en")?`${o.value}-en`:o.value;return[{rel:"canonical",href:`${cn.url}/${p.value}${o.value}`},{rel:"alternate",hreflang:"zh-CN",href:`${cn.url}/zh${h}`},{rel:"alternate",hreflang:"en",href:`${cn.url}/en${d}`},{rel:"alternate",hreflang:"x-default",href:`${cn.url}/zh${h}`}]}),meta:ns(()=>[{property:"og:url",content:`${cn.url}/${p.value}${o.value}`}])});const c=i=>{const h=Lo(i);h&&t.setLocale(h)},u=i=>i instanceof HTMLElement?i.tagName==="INPUT"||i.tagName==="TEXTAREA"||i.isContentEditable:!1,r=i=>{i.key==="/"&&!u(i.target)&&(i.preventDefault(),n.value=!0)};return Vs(()=>e.fullPath,c),En(()=>{c(e.fullPath),window.addEventListener("keydown",r)}),zn(()=>{window.removeEventListener("keydown",r)}),(i,h)=>{const d=Za("router-view");return I(),N("div",yy,[bs(ty,{onOpenSearch:h[0]||(h[0]=b=>n.value=!0)}),S("main",vy,[bs(d,null,{default:jn(({Component:b})=>[bs(_d,{name:"view",mode:"out-in"},{default:jn(()=>[(I(),Aa(Ch(b)))]),_:2},1024)]),_:1})]),bs(gy),bs(_y),n.value?(I(),Aa(ss(a),{key:0,onClose:h[1]||(h[1]=b=>n.value=!1)})):ms("",!0)])}}}),Cy=hn(wy,[["__scopeId","data-v-e6413c2e"]]);function du(s){return Math.max(1,Math.round(Number(s)/300))}const ky={class:"post-item__index num","aria-hidden":"true"},xy={class:"post-item__meta"},Sy={class:"post-item__cat"},Py={class:"post-item__date"},Ey={class:"post-item__words"},Ty={class:"post-item__reading"},Ay={class:"post-item__title"},Ry={class:"post-item__preview"},Ly={class:"post-item__tags"},My=["onClick"],Dy=["aria-label"],Oy={class:"pagination__side"},Iy=["aria-label"],Ny={class:"pagination__pages"},Fy={key:1,class:"page-ellipsis"},$y=["onClick"],By={key:2,class:"page-ellipsis"},qy={class:"pagination__side pagination__side--right"},zy=["aria-label"],Vy=Ns({__name:"PostList",props:{docs:{},perPage:{default:6}},setup(s){const n=s,{t:a,locale:e}=Us(),t=ya(),l=qe(),o=rs(1),p=rs(5),{data:c}=fa(()=>v_(),[]),{data:u}=fa(()=>Gt(),[]),r=ns(()=>a("pagination")),i=ns(()=>Math.max(1,Math.ceil(n.docs.length/n.perPage))),h=ns(()=>{const z=(o.value-1)*n.perPage,Q=z+n.perPage;return n.docs.slice(z,Q)}),d=ns(()=>{const z=[],Q=i.value,cs=o.value,ls=p.value;if(Q<=ls){for(let us=1;us<=Q;us++)z.push(us);return z}z.push(cs);let X=1;for(;z.length<ls&&(cs-X>=1&&!z.includes(cs-X)&&z.push(cs-X),!(z.length>=ls));)cs+X<=Q&&!z.includes(cs+X)&&z.push(cs+X),X++;return z.includes(1)||(z.pop(),z.push(1)),z.length<ls&&!z.includes(Q)?z.push(Q):z.length>=ls&&!z.includes(Q)&&(z.pop(),z.push(Q)),z.sort((us,Fs)=>us-Fs)}),b=ns(()=>d.value.filter(z=>z!==1&&z!==i.value)),f=ns(()=>d.value.includes(1)),C=ns(()=>d.value.includes(i.value)),m=ns(()=>i.value>p.value&&d.value[1]>2),g=ns(()=>i.value>p.value&&d.value[d.value.length-2]<i.value-1),j=ns(()=>{const z={};return u.value.forEach(Q=>Q.items.forEach(cs=>cs.categories.forEach(ls=>{ls.key&&ls.title&&(z[ls.key]=ls.title)}))),z});function y(z){const Q=c.value.find(ls=>ls.title===z.title);let cs;if(Q&&Q.relativePath)cs=`notes/${Q.relativePath}.md`;else{const ls=z.category[1]||"notes",X=z.title.toLowerCase().replace(/[^a-z0-9]/g,"-");cs=`${ls}/${X}.md`,e.value==="en-US"&&(cs=cs.replace(".md","-en.md"))}return Jn(`/article/${cs.replace(/\.md$/,"")}`)}function x(z){o.value=z;const Q={...t.query,page:String(z)};l.push({path:t.path,query:Q}).catch(()=>{}),Qa()}function R(z){z>=1&&z<=i.value&&x(z)}function E(){o.value>1&&x(o.value-1)}function O(){o.value<i.value&&x(o.value+1)}function P(z){Oo(l,z,{locale:ea(e.value),query:{...t.query,page:"1"}})}function U(){p.value=window.innerWidth<480?3:5}function Y(z){const[Q,cs]=z;if(!cs)return Q;const ls=j.value[cs]||cs;return`${Q} / ${ls}`}function V(z){return a("postReadingTime",{minutes:du(z)})}return Vs(()=>n.docs,()=>{o.value=1}),Vs(()=>t.query.page,z=>{const Q=parseInt(String(z)),cs=Number.isFinite(Q)&&Q>=1?Math.min(Q,i.value):1;cs!==o.value&&(o.value=cs,Qa())}),En(()=>{const z=parseInt(String(t.query.page));o.value=Number.isFinite(z)&&z>=1?Math.min(z,i.value):1,U(),window.addEventListener("resize",U)}),zn(()=>{window.removeEventListener("resize",U)}),(z,Q)=>{const cs=Za("router-link"),ls=Ne("reveal");return I(),N("div",null,[(I(!0),N(gs,null,As(h.value,(X,us)=>an((I(),N("article",{key:X.id,class:"post-item",style:qn({"--reveal-delay":Math.min(us,5)*40+"ms"})},[S("span",ky,G(String((o.value-1)*n.perPage+us+1).padStart(2,"0")),1),bs(cs,{to:y(X),class:"post-item__main"},{default:jn(()=>[S("div",xy,[S("span",Sy,G(Y(X.category)),1),Q[5]||(Q[5]=S("span",{class:"divider-v"},null,-1)),S("span",Py,[Q[2]||(Q[2]=S("i",{class:"fas fa-calendar-alt"},null,-1)),Cs(G(X.date),1)]),Q[6]||(Q[6]=S("span",{class:"divider-v"},null,-1)),S("span",Ey,[Q[3]||(Q[3]=S("i",{class:"fas fa-file-lines"},null,-1)),Cs(G(X.wordCount)+" "+G(ss(a)("wordUnit")),1)]),Q[7]||(Q[7]=S("span",{class:"divider-v"},null,-1)),S("span",Ty,[Q[4]||(Q[4]=S("i",{class:"fas fa-clock"},null,-1)),Cs(G(V(X.wordCount)),1)])]),S("h3",Ay,G(X.title),1),S("p",Ry,G(X.preview),1)]),_:2},1032,["to"]),S("div",Ly,[(I(!0),N(gs,null,As(X.tags,Fs=>(I(),N("span",{key:Fs,class:"post-item__tag",onClick:Qs=>P(Fs)},G(Fs),9,My))),128))]),bs(cs,{to:y(X),class:"post-item__arrow","aria-label":X.title},{default:jn(()=>[...Q[8]||(Q[8]=[S("i",{class:"fas fa-arrow-right"},null,-1)])]),_:1},8,["to","aria-label"])],4)),[[ls]])),128)),i.value>1?(I(),N("nav",{key:0,class:"pagination","aria-label":r.value},[S("div",Oy,[o.value>1?(I(),N("button",{key:0,class:"page-btn page-btn--nav",onClick:E,"aria-label":ss(a)("prevPage")},[...Q[9]||(Q[9]=[S("i",{class:"fas fa-chevron-left"},null,-1)])],8,Iy)):ms("",!0)]),S("div",Ny,[f.value?(I(),N("button",{key:0,class:Ts(["page-btn",{"page-btn--active":o.value===1}]),onClick:Q[0]||(Q[0]=X=>R(1))}," 1 ",2)):ms("",!0),m.value?(I(),N("span",Fy,"...")):ms("",!0),(I(!0),N(gs,null,As(b.value,X=>(I(),N("button",{key:X,class:Ts(["page-btn",{"page-btn--active":o.value===X}]),onClick:us=>R(X)},G(X),11,$y))),128)),g.value?(I(),N("span",By,"...")):ms("",!0),C.value&&i.value>1?(I(),N("button",{key:3,class:Ts(["page-btn",{"page-btn--active":o.value===i.value}]),onClick:Q[1]||(Q[1]=X=>R(i.value))},G(i.value),3)):ms("",!0)]),S("div",qy,[o.value<i.value?(I(),N("button",{key:0,class:"page-btn page-btn--nav",onClick:O,"aria-label":ss(a)("nextPage")},[...Q[10]||(Q[10]=[S("i",{class:"fas fa-chevron-right"},null,-1)])],8,zy)):ms("",!0)])],8,Dy)):ms("",!0)])}}}),Uy=hn(Vy,[["__scopeId","data-v-73af845d"]]),Hy={class:"page-section home-section"},Gy={class:"hero"},Wy={class:"hero__main"},Ky={class:"hero__greeting"},Xy={class:"hero__greeting-mark"},Yy={class:"hero__name"},Qy={class:"hero__bio"},Jy={class:"hero__stats"},Zy={class:"hero__stat"},sv={class:"hero__stat-num num"},nv={class:"hero__stat-label"},av={class:"hero__stat"},ev={class:"hero__stat-num num"},tv={class:"hero__stat-label"},lv={class:"hero__stat"},ov={class:"hero__stat-num num"},pv={class:"hero__stat-label"},cv={key:0,class:"tags-row"},rv=["onClick"],iv={class:"posts-header"},uv={key:0,class:"posts-header__tag-title"},hv={class:"posts-header__actions"},dv=["aria-label"],mv=["aria-label"],jv=Ns({name:"HomeView",__name:"Home",setup(s){It({title:"zorrooz’s blog - Home"});const{t:n,locale:a}=Us(),e=ya(),t=qe(),{data:l}=fa(()=>y_(),[]),{data:o}=fa(()=>w_(),[]),p=cn.author,c=ns(()=>typeof e.query.tag=="string"?e.query.tag:""),u=ns(()=>typeof e.query.from=="string"?e.query.from:""),r=ns(()=>{const j=c.value;return j?l.value.filter(y=>y.tags.includes(j)):l.value}),i=ns(()=>{const j=new Map;l.value.forEach(O=>O.tags.forEach(P=>j.set(P,(j.get(P)||0)+1)));const y=Array.from(j.entries()).map(([O,P])=>({name:O,count:P})).sort((O,P)=>O.name.localeCompare(P.name)),x=y.length,R=Math.ceil(x/3),E=new Map([...y].sort((O,P)=>P.count-O.count).map((O,P)=>[O.name,P]));return y.map(O=>{const P=E.get(O.name)??0;return{...O,level:P<R?"lg":P<R*2?"md":"sm"}})}),h=ns(()=>l.value.length),d=ns(()=>o.value.length),b=ns(()=>l.value.reduce((j,y)=>j+y.wordCount,0)),f=ns(()=>{const j=b.value;return j>=1e6?(j/1e6).toFixed(j%1e6?1:0)+"M":j>=1e3?(j/1e3).toFixed(j%1e3?1:0)+"K":String(j)});function C(){const j={...e.query};delete j.tag,delete j.from,j.page="1",t.push({path:e.path,query:j}).catch(()=>{}),Qa()}function m(){u.value&&t.push(u.value).catch(()=>{})}function g(j){Oo(t,j,{locale:ea(a.value),query:{...e.query,page:"1"}})}return Vs(a,(j,y)=>{j!==y&&c.value&&C()}),(j,y)=>{const x=Ne("reveal");return I(),N("div",Hy,[S("header",Gy,[S("div",Wy,[S("span",Ky,[S("span",Xy,G(ss(n)("greetingPrefix")),1),Cs(G(ss(n)("greeting")),1)]),S("h1",Yy,G(ss(p)),1),S("p",Qy,G(ss(n)("developer")),1)]),S("div",Jy,[S("div",Zy,[S("span",sv,G(h.value),1),S("span",nv,G(ss(n)("articles")),1)]),S("div",av,[S("span",ev,G(d.value),1),S("span",tv,G(ss(n)("tags")),1)]),S("div",lv,[S("span",ov,G(f.value),1),S("span",pv,G(ss(n)("words")),1)])])]),c.value?ms("",!0):an((I(),N("div",cv,[(I(!0),N(gs,null,As(i.value,R=>(I(),N("span",{key:R.name,class:Ts(["tag",`tag--${R.level}`]),onClick:E=>g(R.name)},G(R.name),11,rv))),128))])),[[x]]),an((I(),N("div",iv,[S("h2",{class:Ts(["posts-header__title",{"posts-header__title--tag":c.value}])},[c.value?(I(),N("span",uv,"# "+G(c.value),1)):(I(),N(gs,{key:1},[Cs(G(ss(n)("recentPosts")),1)],64))],2),S("div",hv,[u.value?(I(),N("button",{key:0,class:"chip-close",onClick:m,"aria-label":ss(n)("backToArticle")},[y[0]||(y[0]=S("i",{class:"fas fa-arrow-left"},null,-1)),Cs(G(ss(n)("backToArticle")),1)],8,dv)):ms("",!0),c.value?(I(),N("button",{key:1,class:"chip-close",onClick:C,"aria-label":ss(n)("close")},[y[1]||(y[1]=S("i",{class:"fas fa-times"},null,-1)),Cs(G(ss(n)("clearFilter")),1)],8,mv)):ms("",!0)])])),[[x]]),bs(Uy,{docs:r.value,perPage:5},null,8,["docs"])])}}}),gv=hn(jv,[["__scopeId","data-v-774fc76a"]]),fv={class:"page-section category-view"},bv={class:"category-head"},_v={class:"article-title"},yv={class:"cat-section__header"},vv={class:"cat-section__title"},wv={class:"cat-section__count"},Cv={class:"cat-grid"},kv={class:"cat-card__head"},xv={class:"cat-card__name"},Sv={class:"cat-card__ext-links"},Pv=["href"],Ev=["href"],Tv={class:"cat-card__desc"},Av={key:0,class:"cat-card__stats"},Rv={key:0,class:"cat-stat"},Lv={key:1,class:"cat-stat"},Mv={key:2,class:"cat-stat"},Dv={key:1,class:"cat-card__stats"},Ov={key:0,class:"cat-stat"},Iv={key:1,class:"cat-stat"},Nv={key:2,class:"cat-stat"},Fv={key:2,class:"cat-card__tags"},$v={class:"cat-card__links"},Bv=["onClick"],qv=Ns({name:"CategoryView",__name:"Category",setup(s){It({title:"zorrooz’s blog - Categories"});const{t:n}=Us(),a=qe(),{data:e}=fa(()=>Gt(),[]),t=ns(()=>n("categories")),l=ns(()=>n("notes")),o=ns(()=>n("projects")),p=ns(()=>n("topics")),c=ns(()=>n("seeMore"));function u(C){return C===l.value?"fa-book-open":C===o.value?"fa-folder-open":C===p.value?"fa-flask":"fa-folder"}function r(C){const m=C.items;if(C.title===l.value){const g=m.reduce((j,y)=>j+y.stats.postsCount,0);return n("countPosts",{count:g})}return C.title===o.value?n("countProjects",{count:m.length}):n("countTopics",{count:m.length})}function i(C){return C.stats.latestDate||""}function h(C){return!C||!C.trim()?"":/^https?:\/\//i.test(C)?C:"https://"+C.replace(/^\/+/,"")}function d(C){return!C||!C.trim()?"":/^https?:\/\//i.test(C)?C:/^10\.\d{4,9}\//.test(C)?"https://doi.org/"+C:"https://"+C.replace(/^\/+/,"")}function b(C){return!!(C.url||C.github||C.doi)}function f(C){if(C.root){a.push(Jn(C.root)).catch(g=>{g.name!=="NavigationDuplicated"&&!g.toString().includes("Navigation cancelled")&&console.error("Navigation error:",g)});return}const m=["url","github","doi"];for(const g of m){const j=C[g];if(!j)continue;const y=g==="doi"?d(j):h(j);if(y){window.open(y,"_blank","noopener,noreferrer");return}}}return(C,m)=>{const g=Ne("reveal");return I(),N("div",fv,[S("div",bv,[an((I(),N("h1",_v,[Cs(G(t.value),1)])),[[g]])]),(I(!0),N(gs,null,As(ss(e),(j,y)=>(I(),N("div",{key:y,class:"cat-section"},[an((I(),N("div",yv,[S("h2",vv,[S("i",{class:Ts([["fas",u(j.title)],"cat-section__icon"]),"aria-hidden":"true"},null,2),Cs(" "+G(j.title),1)]),S("span",wv,G(r(j)),1)])),[[g]]),S("div",Cv,[(I(!0),N(gs,null,As(j.items,(x,R)=>an((I(),N("div",{key:R,class:"cat-card",style:qn({"--reveal-delay":Math.min(Number(R),5)*40+"ms"})},[S("div",kv,[S("div",xv,G(x.title),1),S("div",Sv,[x.github?(I(),N("a",{key:0,href:h(x.github),target:"_blank",rel:"noopener noreferrer",class:"cat-card__ext-link","aria-label":"GitHub"},[...m[0]||(m[0]=[S("i",{class:"fab fa-github"},null,-1)])],8,Pv)):ms("",!0),x.doi?(I(),N("a",{key:1,href:d(x.doi),target:"_blank",rel:"noopener noreferrer",class:"cat-card__ext-link","aria-label":"DOI"},[...m[1]||(m[1]=[S("i",{class:"fas fa-link"},null,-1)])],8,Ev)):ms("",!0)])]),S("p",Tv,G(x.desc),1),j.title===l.value?(I(),N("div",Av,[x.stats?.postsCount?(I(),N("span",Rv,[m[2]||(m[2]=S("i",{class:"fas fa-file-lines"},null,-1)),Cs(G(ss(n)("countPosts",{count:x.stats.postsCount})),1)])):ms("",!0),x.stats?.totalWords?(I(),N("span",Lv,[m[3]||(m[3]=S("i",{class:"fas fa-font"},null,-1)),Cs(G(ss(n)("countWords",{count:x.stats.totalWords})),1)])):ms("",!0),i(x)?(I(),N("span",Mv,[m[4]||(m[4]=S("i",{class:"fas fa-clock"},null,-1)),Cs(G(i(x)),1)])):ms("",!0)])):(I(),N("div",Dv,[x.language?(I(),N("span",Ov,[m[5]||(m[5]=S("i",{class:"fas fa-code"},null,-1)),Cs(G(x.language),1)])):ms("",!0),x.year?(I(),N("span",Iv,[m[6]||(m[6]=S("i",{class:"fas fa-calendar"},null,-1)),Cs(G(x.year),1)])):ms("",!0),x.license?(I(),N("span",Nv,[m[7]||(m[7]=S("i",{class:"fas fa-scale-balanced"},null,-1)),Cs(G(x.license),1)])):ms("",!0)])),Array.isArray(x.tags)&&x.tags.length?(I(),N("div",Fv,[(I(!0),N(gs,null,As(x.tags,(E,O)=>(I(),N("span",{key:O,class:"cat-card__tag"},G(E),1))),128))])):ms("",!0),S("div",$v,[b(x)||x.root?(I(),N("a",{key:0,class:"cat-card__link",onClick:vn(E=>f(x),["prevent"])},[Cs(G(c.value),1),m[8]||(m[8]=S("i",{class:"fas fa-arrow-right"},null,-1))],8,Bv)):ms("",!0)])],4)),[[g]])),128))])]))),128))])}}}),zv=hn(qv,[["__scopeId","data-v-b3912557"]]),Vv={class:"page-section resource-view"},Uv={class:"resource-head"},Hv={class:"article-title"},Gv={class:"resource-subtitle"},Wv={class:"res-layout"},Kv={class:"res-sidebar"},Xv={class:"res-group__label"},Yv={class:"res-group__count num"},Qv={class:"res-group__items"},Jv=["onClick"],Zv={class:"res-main"},sw={class:"res-groups"},nw={class:"res-group-block__title"},aw={class:"res-grid"},ew=["href"],tw={class:"res-card__head"},lw={class:"res-card__name"},ow={class:"res-card__ext-links"},pw={key:0,class:"res-card__ext-link","aria-label":"DOI"},cw={class:"res-card__desc"},rw={class:"res-card__footer"},iw={class:"res-card__url"},uw=Ns({name:"ResourceView",__name:"Resource",setup(s){const{t:n}=Us(),{data:a}=fa(()=>C_(),[]),e=rs(null),t=ns(()=>n("resources")),l=ns(()=>n("resourceSubtitle")),o=ns(()=>{const i=e.value;return i?i.children?.length?i.children.filter(h=>h.items?.length):i.items?.length?[{title:i.title,items:i.items}]:[]:[]});function p(i){e.value=i}function c(i){return e.value===i}function u(i){return i?i.replace(/^https?:\/\//,"").replace(/\/$/,""):""}function r(i){return!!i&&(i.includes("doi.org")||/^10\.\d{4,9}\//.test(i))}return Vs(a,i=>{e.value=i[0]?.children?.[0]||null},{immediate:!0}),(i,h)=>{const d=Ne("reveal");return I(),N("div",Vv,[S("header",Uv,[an((I(),N("h1",Hv,[Cs(G(t.value),1)])),[[d]]),S("p",Gv,[h[0]||(h[0]=S("i",{class:"fas fa-circle-info resource-head__icon"},null,-1)),Cs(G(l.value),1)])]),S("div",Wv,[S("aside",Kv,[(I(!0),N(gs,null,As(ss(a),b=>(I(),N("div",{key:b.title,class:"res-group"},[S("div",Xv,[S("span",null,G(b.title),1),S("span",Yv,G(b.children?.length||0),1)]),S("div",Qv,[(I(!0),N(gs,null,As(b.children,f=>(I(),N("button",{key:f.title,class:Ts(["res-item",{"res-item--active":c(f)}]),onClick:C=>p(f)},G(f.title),11,Jv))),128))])]))),128))]),h[3]||(h[3]=S("div",{class:"res-divider","aria-hidden":"true"},null,-1)),S("main",Zv,[S("div",sw,[(I(!0),N(gs,null,As(o.value,b=>(I(),N("section",{key:b.title,class:"res-group-block"},[S("h3",nw,G(b.title),1),S("div",aw,[(I(!0),N(gs,null,As(b.items,(f,C)=>an((I(),N("a",{key:f.name,href:f.url,target:"_blank",rel:"noopener noreferrer",class:"res-card",style:qn({"--reveal-delay":Math.min(Number(C),5)*40+"ms"})},[S("div",tw,[S("span",lw,G(f.name),1),S("span",ow,[r(f.url)?(I(),N("span",pw,[...h[1]||(h[1]=[S("i",{class:"fas fa-link"},null,-1)])])):ms("",!0)])]),S("p",cw,G(f.desc),1),S("div",rw,[S("span",iw,G(u(f.url)),1),h[2]||(h[2]=S("i",{class:"fas fa-arrow-up-right-from-square res-card__arrow"},null,-1))])],12,ew)),[[d]])),128))])]))),128))])])])])}}}),hw=hn(uw,[["__scopeId","data-v-30eaadfe"]]),dw="/assets/avatar-DQvqWlfS.png",mw={class:"page-section about-view"},jw={class:"about-head"},gw={class:"about-head__identity"},fw={class:"about-head__avatar"},bw=["src"],_w={key:1,class:"about-head__initial"},yw={class:"about-head__names"},vw={class:"about-head__name"},ww={class:"about-head__role"},Cw={key:0,class:"about-intro"},kw={class:"about-body"},xw={key:0,class:"about-section"},Sw={class:"about-section__title"},Pw={class:"timeline"},Ew={class:"tl-year num"},Tw={class:"tl-content"},Aw={class:"tl-title"},Rw={key:0,class:"tl-desc"},Lw={key:1,class:"about-grid"},Mw={class:"about-cell__title"},Dw={key:0,class:"about-cell__item-name"},Ow={key:1,class:"about-cell__item-desc"},Iw={class:"about-foot"},Nw={class:"about-foot__contacts"},Fw=["href"],$w={class:"about-foot__thanks"},Bw=Ns({name:"AboutView",__name:"About",setup(s){const{t:n}=Us(),{data:a}=fa(()=>k_(),iu),e=Object.assign({"/src/assets/avatar.png":dw}),t=ns(()=>Object.values(e)[0]||""),l=ns(()=>n("thanks")),o=cn.author,p=o.trim().charAt(0).toUpperCase(),c=ns(()=>a.value.introduction),u=ns(()=>a.value.experience),r=ns(()=>a.value.section);return(i,h)=>{const d=Ne("reveal");return I(),N("div",mw,[S("header",jw,[an((I(),N("div",gw,[S("div",fw,[t.value?(I(),N("img",{key:0,src:t.value,alt:"avatar",class:"about-head__avatar-img",draggable:"false"},null,8,bw)):(I(),N("span",_w,G(ss(p)),1))]),S("div",yw,[S("h1",vw,G(ss(o)),1),S("p",ww,G(ss(n)("developer")),1)])])),[[d]]),c.value?an((I(),N("p",Cw,[Cs(G(c.value),1)])),[[d]]):ms("",!0)]),S("main",kw,[u.value.length?an((I(),N("section",xw,[S("div",Sw,G(ss(n)("experience")),1),S("div",Pw,[(I(!0),N(gs,null,As(u.value,(b,f)=>(I(),N("div",{key:f,class:"tl-item"},[S("div",Ew,G(b.year),1),S("div",Tw,[S("div",Aw,G(b.title),1),b.desc?(I(),N("div",Rw,G(b.desc),1)):ms("",!0)])]))),128))])])),[[d]]):ms("",!0),r.value.length?(I(),N("div",Lw,[(I(!0),N(gs,null,As(r.value,(b,f)=>an((I(),N("div",{key:f,class:"about-cell",style:qn({"--reveal-delay":Math.min(Number(f),5)*40+"ms"})},[S("div",Mw,G(b.title),1),(I(!0),N(gs,null,As(b.items,(C,m)=>(I(),N("div",{key:m,class:"about-cell__item"},[C.name?(I(),N("span",Dw,G(C.name),1)):ms("",!0),C.desc?(I(),N("span",Ow,G(C.desc),1)):ms("",!0)]))),128))],4)),[[d]])),128))])):ms("",!0)]),S("footer",Iw,[S("div",Nw,[(I(!0),N(gs,null,As(ss(a).contacts,b=>(I(),N("a",{key:b.label,href:b.link,target:"_blank",rel:"noopener noreferrer",class:"about-foot__contact"},[S("i",{class:Ts(b.icon)},null,2),S("span",null,G(b.value),1)],8,Fw))),128))]),S("p",$w,G(l.value),1)])])}}}),qw=hn(Bw,[["__scopeId","data-v-e9b50f8b"]]);async function mu(s){try{await navigator.clipboard.writeText(s)}catch{const n=document.createElement("textarea");n.value=s,document.body.appendChild(n),n.select(),document.execCommand("copy"),document.body.removeChild(n)}}const zw=["innerHTML"],Mc=`
  <svg width="16" height="16" viewBox="0 0 14 14" fill="currentColor">
    <path d="M3 2C2.44772 2 2 2.44772 2 3V9C2 9.55228 2.44772 10 3 10H9C9.55228 10 10 9.55228 10 9V3C10 2.44772 9.55228 2 9 2H3ZM1 3C1 1.89543 1.89543 1 3 1H9C10.1046 1 11 1.89543 11 3V9C11 10.1046 10.1046 11 9 11H3C1.89543 11 1 10.1046 1 9V3Z"/>
    <path d="M5 4C4.44772 4 4 4.44772 4 5V11C4 11.5523 4.44772 12 5 12H11C11.5523 12 12 11.5523 12 11V5C12 4.44772 11.5523 4 11 4H5Z"/>
  </svg>
`,Vw=`
  <svg width="16" height="16" viewBox="0 0 14 14" fill="currentColor">
    <path d="M11.3536 3.64645C11.5488 3.84171 11.5488 4.15829 11.3536 4.35355L5.35355 10.3536C5.15829 10.5488 4.84171 10.5488 4.64645 10.3536L2.64645 8.35355C2.45118 8.15829 2.45118 7.84171 2.64645 7.64645C2.84171 7.45118 3.15829 7.45118 3.35355 7.64645L5 9.29289L10.6464 3.64645C10.8417 3.45118 11.1583 3.45118 11.3536 3.64645Z"/>
  </svg>
`,Uw=Ns({__name:"RenderMarkdown",props:{rawMarkdown:{default:""},articlePath:{default:""},articleTitle:{default:""}},emits:["markdown-rendered"],setup(s,{emit:n}){const{t:a}=Us(),e=Object.assign({}),t=s,l=n,o=rs(""),p=at("markdownContainer");async function c(g){const j=u(g,t.articlePath);o.value=j,await bn(),l("markdown-rendered"),d(),b(),r()}function u(g,j){try{const y=j.replace(/^[./]*/,"").replace(/\.md$/,"").split("/").slice(0,-1).join("/"),x=R=>{if(/^(https?:)?\/\//i.test(R)||R.startsWith("/"))return R;const E=(y+"/"+R).split("/").filter(Y=>Y&&Y!=="."),O=[];E.forEach(Y=>Y===".."?O.pop():O.push(Y));const P=O.join("/"),U=[`@data/content-src/${P}`,`${P}`];for(const Y of U){const V=Object.keys(e).find(z=>z.endsWith(`/${P}`)||z===Y);if(V)return e[V]}return R};return g.replace(/<img\s+([^>]*?)src=["']([^"']+)["'](.*?)>/gi,(R,E,O,P)=>`<img ${E}src="${x(O.trim())}"${P}>`)}catch(y){return console.warn("rewriteImageLinks failed",y),g}}function r(){const g=p.value;g&&(i(g),h(g))}function i(g){if(!t.articleTitle)return;const j=t.articleTitle.trim().toLowerCase();g.querySelectorAll("h1").forEach(y=>{y.textContent.trim().toLowerCase()===j?y.remove():y.replaceWith(Object.assign(document.createElement("h2"),{...Object.fromEntries(Array.from(y.attributes).map(R=>[R.name,R.value])),innerHTML:y.innerHTML}))})}function h(g){const j=(y,x)=>{const R=window.scrollY+y.getBoundingClientRect().top-88;window.scrollTo({top:Math.max(0,R),behavior:"smooth"}),setTimeout(()=>x.blur(),300)};g.querySelectorAll("h2, h3, h4, h5, h6").forEach(y=>{y.querySelector(".heading-anchor")?.remove();const x=Object.assign(document.createElement("button"),{type:"button",className:"heading-anchor",textContent:"#",ariaLabel:a("anchorHeading"),tabIndex:0,ariaHidden:"false"});x.addEventListener("click",R=>{R.stopPropagation(),j(y,x)}),x.addEventListener("keydown",R=>{(R.key==="Enter"||R.key===" ")&&(R.preventDefault(),j(y,x))}),y.appendChild(x)})}function d(){const g=p.value;g&&g.querySelectorAll("pre").forEach(j=>{if(j.querySelector(".code-block-header"))return;const y=j.querySelector("code");if(!y)return;const x=(y.className.match(/language-(\w+)/)||["","text"])[1],R=document.createElement("div");R.className="code-block-header d-flex align-items-center justify-content-between";const E=document.createElement("span");E.className="code-language",E.textContent=x;const O=document.createElement("button");O.type="button",O.className="copy-button btn-icon d-flex align-items-center justify-content-center",O.setAttribute("aria-label",a("copyCode")),O.innerHTML=Mc,O.addEventListener("click",()=>C(y.textContent??"",O)),R.append(E,O);const P=document.createElement("div");P.className="code-block-wrapper",j.parentNode?.insertBefore(P,j),P.append(R,j)})}function b(){const g=p.value;g&&g.querySelectorAll("table").forEach(j=>{if(j.closest(".table-copyable"))return;const y=document.createElement("div");y.className="table-copyable";const x=document.createElement("button");x.type="button",x.className="table-copy-btn",x.setAttribute("aria-label",a("copyTable")),x.innerHTML=Mc,x.addEventListener("click",()=>f(j,x)),y.append(x),j.parentNode?.insertBefore(y,j),y.append(j)})}function f(g,j){const x=Array.from(g.querySelectorAll("tr")).map(R=>Array.from(R.querySelectorAll("th, td")).map(E=>(E.textContent||"").trim().replace(/\s+/g," ")).join("	"));C(x.join(`
`),j)}async function C(g,j){try{await mu(g)}catch(y){console.error(a("copyFailed"),y)}finally{m(j)}}function m(g){const j=g.innerHTML;g.style.color="var(--primary)",g.innerHTML=Vw,setTimeout(()=>{g.innerHTML=j,g.style.color=""},1200)}return Vs(()=>t.rawMarkdown,g=>c(g),{immediate:!0}),(g,j)=>(I(),N("div",{class:"markdown-body",innerHTML:o.value,ref_key:"markdownContainer",ref:p},null,8,zw))}}),Hw={class:"on-this-page"},Gw={class:"otp-header"},Ww={class:"otp-title"},Kw={class:"otp-list"},Xw=["onClick","onKeydown"],Yw={class:"otp-text"},Qw={key:0,class:"otp-sublist"},Jw=["onClick","onKeydown"],Zw={class:"otp-subtext"},s1=36,n1=Ns({__name:"OnThisPage",props:{containerSelector:{default:".markdown-body"},levels:{default:()=>[2,3,4,5,6]},offset:{default:8}},emits:["navigate"],setup(s,{expose:n,emit:a}){const e=s,t=a,{t:l}=Us(),o=rs([]),p=rs(""),c=rs(null),u=rs(null),r=rs(null),i=rs(null),h=rs(null),d=rs(null);function b(){r.value&&(r.value.disconnect(),r.value=null),i.value&&(clearTimeout(i.value),i.value=null),h.value&&(clearInterval(h.value),h.value=null),d.value&&(d.value.disconnect(),d.value=null)}function f(){o.value=[],p.value="",b(),bn(()=>g())}function C(){y(),x(),bn(m)}function m(){if(!u.value||!c.value)return;const E=c.value.querySelector(".otp-list");if(!E)return;const O=E.querySelector(".otp-item.active > .otp-link, .otp-subitem.active > .otp-sublink");if(!O){u.value.style.opacity="0";return}u.value.style.top=`${O.offsetTop+s1}px`,u.value.style.opacity="1"}Vs(p,()=>bn(m));function g(){setTimeout(()=>C(),100);const E=()=>{const O=document.querySelector(e.containerSelector);O&&(h.value&&(clearInterval(h.value),h.value=null),r.value||(r.value=new MutationObserver(()=>{clearTimeout(i.value??void 0),i.value=window.setTimeout(()=>C(),100)}),r.value.observe(O,{childList:!0,subtree:!0,attributes:!0,attributeFilter:["id"]})),C())};E(),!h.value&&!document.querySelector(e.containerSelector)&&(h.value=window.setInterval(E,200))}function j(E){try{const O=E.cloneNode(!0);return O.querySelectorAll(".heading-anchor")?.forEach(P=>P.remove()),(O.textContent||"").replace(/\s*#\s*$/,"").trim()}catch{return(E.textContent||"").replace(/\s*#\s*$/,"").trim()}}function y(){const E=document.querySelector(e.containerSelector);if(!E){o.value=[];return}const O=e.levels.map(cs=>`h${cs}`).join(","),P=Array.from(E.querySelectorAll(O));if(P.length===0){o.value=[];return}P.forEach(cs=>{if(!cs.id){let ls=cs.textContent.trim().toLowerCase().replace(/[^\u4e00-\u9fa5a-zA-Z0-9\s-]/g,"").replace(/\s+/g,"-");ls||(ls=`section-${Math.random().toString(36).substring(2,9)}`);let X=ls,us=1;for(;document.getElementById(X);)X=`${ls}-${us++}`;cs.id=X}});const U=new Set(e.levels),Y=Math.min(...e.levels),V=Y+1,z=[];let Q=null;for(const cs of P){const ls=parseInt(cs.tagName.substring(1),10);if(!U.has(ls))continue;const X={id:cs.id,text:j(cs),level:ls,children:[]};ls===Y?(z.push(X),Q=X):Q&&ls>=V?Q.children.push(X):z.push(X)}o.value=z}function x(){d.value&&(d.value.disconnect(),d.value=null);const E=document.querySelector(e.containerSelector);if(!E)return;const O=e.levels.map(U=>`h${U}`).join(","),P=Array.from(E.querySelectorAll(O));P.length!==0&&(d.value=new IntersectionObserver(U=>{for(const Y of U)Y.isIntersecting&&(p.value=Y.target.id)},{rootMargin:`-${e.offset}px 0px -60% 0px`,threshold:0}),P.forEach(U=>d.value?.observe(U)))}function R(E){t("navigate",E);const O=document.getElementById(E);if(!O)return;const P=O.getBoundingClientRect().top+window.scrollY-e.offset,U=typeof window.matchMedia=="function"&&window.matchMedia("(prefers-reduced-motion: reduce)").matches,Y=()=>{window.scrollTo({top:P,behavior:U?"auto":"smooth"})};try{(document.body&&document.body.style&&document.body.style.overflow)==="hidden"?setTimeout(Y,80):Y()}catch{Y()}}return En(()=>{y(),x(),bn(()=>{m(),g()})}),zn(()=>{b()}),n({refreshToc:C,resetToc:f}),(E,O)=>(I(),N("nav",Hw,[S("div",{class:"otp-content",ref_key:"otpContent",ref:c},[S("span",{class:"otp-marker",ref_key:"marker",ref:u,"aria-hidden":"true"},null,512),S("div",Gw,[S("span",Ww,G(ss(l)("tableOfContents")),1)]),S("ul",Kw,[(I(!0),N(gs,null,As(o.value,(P,U)=>(I(),N("li",{key:U,class:Ts(["otp-item",{active:p.value===P.id}])},[S("a",{class:"otp-link",role:"button",tabindex:"0",onClick:vn(Y=>R(P.id),["prevent"]),onKeydown:Sp(vn(Y=>R(P.id),["prevent"]),["enter"])},[S("span",Yw,G(P.text),1)],40,Xw),P.children&&P.children.length?(I(),N("ul",Qw,[(I(!0),N(gs,null,As(P.children,(Y,V)=>(I(),N("li",{key:V,class:Ts(["otp-subitem",{active:p.value===Y.id}])},[S("a",{class:"otp-sublink",role:"button",tabindex:"0",onClick:vn(z=>R(Y.id),["prevent"]),onKeydown:Sp(vn(z=>R(Y.id),["prevent"]),["enter"])},[S("span",Zw,G(Y.text),1)],40,Jw)],2))),128))])):ms("",!0)],2))),128))])],512)]))}}),ju=hn(n1,[["__scopeId","data-v-02208124"]]),a1=["aria-label"],e1={class:"offcanvas-panel offcanvas-right border-start rounded-0 shadow-sm"},t1={class:"offcanvas-section"},l1={class:"offcanvas-card"},o1=Ns({__name:"TocDrawer",setup(s){const{t:n}=Us(),a=rs(!1),e=rs(!1),t=rs(null),{isDragging:l,buttonTop:o,clampTop:p,dispatchBaseTop:c,onTouchStart:u,onTouchMove:r,onTouchEnd:i,subscribe:h,unsubscribe:d}=hu({sourceId:"toc",defaultTop:(typeof window<"u"?window.innerHeight:1024)-160,mode:"stack",onRelease:f});function b(){const g=window.innerWidth<992;a.value=g,o.value=p(o.value)}function f(){l.value||(t.value=ib(),e.value=!0)}function C(){e.value=!1,Tc(t.value)}function m(){C()}return En(()=>{window.addEventListener("resize",b,{passive:!0}),h(),b(),c()}),zn(()=>{window.removeEventListener("resize",b),d(),Tc(t.value)}),(g,j)=>(I(),N(gs,null,[an(S("button",{class:"toc-drawer-btn d-lg-none d-flex align-items-center justify-content-center",onClick:f,"aria-label":ss(n)("openToc"),onTouchstart:j[0]||(j[0]=vn((...y)=>ss(u)&&ss(u)(...y),["prevent","stop"])),onTouchmove:j[1]||(j[1]=vn((...y)=>ss(r)&&ss(r)(...y),["prevent","stop"])),onTouchend:j[2]||(j[2]=vn((...y)=>ss(i)&&ss(i)(...y),["prevent","stop"])),style:qn({top:ss(o)+"px"})},[...j[3]||(j[3]=[S("i",{class:"fas fa-bookmark"},null,-1)])],44,a1),[[yt,a.value]]),e.value?(I(),N("div",{key:0,class:"mobile-offcanvas d-lg-none",onClick:vn(C,["self"])},[S("div",e1,[S("div",t1,[S("div",l1,[bs(ju,{containerSelector:".markdown-body",levels:[2,3],offset:88,onNavigate:m})])])]),S("div",{class:"offcanvas-backdrop",onClick:C})])):ms("",!0)],64))}}),p1=hn(o1,[["__scopeId","data-v-f5b1ecf6"]]),c1={class:"container view-container article-view"},r1={class:"row py-4 px-0"},i1={class:"col-12 col-lg-2 order-2 order-lg-1 docs-sidebar-col"},u1={key:0,class:"navigation-container mb-0"},h1={class:"col-12 col-lg-8 order-1 order-lg-2 docs-main-col",ref:"mainContent"},d1={class:"article-panel"},m1={class:"article-panel__body"},j1={class:"article-content"},g1={key:0,class:"article-meta"},f1={class:"article-title"},b1={class:"article-meta__row"},_1={key:0,class:"meta-line"},y1={key:1,class:"meta-line"},v1=["aria-label"],w1={key:0,class:"article-meta__tags"},C1=["onClick"],k1={key:2,class:"article-navigation"},x1={key:0,class:"article-nav-spacer","aria-hidden":"true"},S1={class:"nav-details"},P1={class:"nav-label"},E1={class:"nav-title"},T1={class:"nav-details"},A1={class:"nav-label"},R1={class:"nav-title"},L1={class:"col-12 col-lg-2 order-3 docs-toc-col"},M1={key:0,class:"toc-container mt-0"},D1=Ns({name:"ArticleView",head(){return{title:this.currentPost?`${this.currentPost.title} - zorrooz’s blog`:"zorrooz’s blog",meta:this.currentPost?.description?[{name:"description",content:this.currentPost.description}]:[]}},__name:"Article",setup(s){const{t:n,locale:a}=Us(),e=ya(),t=qe(),l=rs(""),o=rs(""),p=rs([]),c=rs([]),u=rs(typeof window<"u"?window.innerWidth:1024),r=rs(0),i=rs(!1);let h=null,d=!1;const b=at("onThisPageRef"),f=at("leftSidebarContent"),C=at("rightSidebarContent"),m=ns(()=>u.value>=992),g=ns(()=>!!P.value?.path.startsWith("notes/")),j=ns(()=>n("updatedAt")),y=ns(()=>n("prevPage")),x=ns(()=>n("nextPage")),R=os=>os.replace(/\.md$/,""),E=(os,L)=>{const W=R(os),H=R(L);return!!(W===H||H.endsWith("-en")&&W===H.replace(/-en$/,"")||!H.endsWith("-en")&&W===H+"-en")},O=os=>p.value.find(L=>E(L.path,os))||null,P=ns(()=>o.value?O(o.value):null),U=ns(()=>{if(!P.value)return[];const[os,L]=R(P.value.path).split("/"),W=[],H=(J,hs)=>{if(!hs.trim())return;const w=R(Xl(hs)),[k,A]=w.split("/");k!==os||A!==L||W.push({title:J,path:`${w}.md`})};for(const J of c.value)for(const hs of J.items)hs.name===L&&(hs.articles?.forEach(w=>H(w.title,w.articleUrl)),hs.categories.forEach(w=>w.articles.forEach(k=>H(k.title,k.articleUrl))));return W}),Y=ns(()=>{const os=P.value;return os?U.value.findIndex(L=>R(L.path)===R(os.path)):-1}),V=ns(()=>{const os=Y.value;return os>0?U.value[os-1]:null}),z=ns(()=>{const os=Y.value,L=U.value.length-1;return os>=0&&os<L?U.value[os+1]:null}),Q=ns(()=>du(P.value?.wordCount??0));function cs(os){return n("articleReadingTime",{minutes:os})}function ls(os){return Jn(nu(os))}function X(os){Oo(t,os,{locale:ea(a.value),query:{from:e.fullPath},scroll:!1})}async function us(){if(!P.value)return;let os="";try{const L=await x_(P.value.path);L&&(os=L.trim())}catch{}if(!os){const L=document.querySelector(".markdown-body");if(!L)return;const W=L.cloneNode(!0);W.querySelectorAll(".heading-anchor, .code-block-header, .copy-button, .table-copy-btn").forEach(J=>J.remove());const H=(W.innerText||"").trim();os=`${P.value.title}

${H}`}try{await mu(os)}catch(L){console.error(n("copyFailed"),L)}finally{i.value=!0,h&&clearTimeout(h),h=window.setTimeout(()=>{i.value=!1},1200)}}const Fs=(os,L,W)=>{const H=Xl(L.articleUrl);os.push({title:L.title,path:H,date:W,tags:L.tags,wordCount:L.wordCount})};function Qs(){const os=[];try{const L=Gt();c.value=L,L.forEach(W=>{W.items.forEach(H=>{const J=H.stats.latestDate||"";H.articles?.forEach(hs=>Fs(os,hs,J)),H.categories.forEach(hs=>{const w=hs.stats.latestDate||J;hs.articles.forEach(k=>Fs(os,k,w))})})}),p.value=os}catch{p.value=[],c.value=[]}}function Ms(){u.value=window.innerWidth,en()}function $s(){d||(d=!0,requestAnimationFrame(()=>{d=!1;const L=document.documentElement.scrollHeight-window.innerHeight;r.value=L>0?Math.min(100,window.scrollY/L*100):0}))}function dn(){l.value="";try{const os=O(ct(e.params.path));if(!os)throw new Error(`Article not found: ${e.params.path}`);o.value=os.path,l.value=__(o.value),bn(()=>{typeof window>"u"||(Qa(),en(),b.value?.refreshToc())})}catch{l.value=`# Article Not Found

The requested article could not be loaded. Please check the URL.`,bn(()=>{typeof window>"u"||(Qa(),b.value?.refreshToc())})}}function en(){if(typeof window>"u")return;const L=document.querySelector("header")?.offsetHeight||60,W=window.innerHeight,H=Math.max(200,W-L-24-24);[f.value,C.value].forEach(hs=>{hs&&(hs.style.maxHeight=`${H}px`,hs.style.overflowY="auto")})}function Mn(){en(),bn(()=>b.value?.refreshToc())}return Qs(),dn(),Vs(a,(os,L)=>{os!==L&&(Qs(),dn())}),Vs(()=>e.params.path,(os,L)=>{const W=ct(L),H=ct(os);W!==H&&(b.value?.resetToc(),dn())}),Vs(l,()=>{bn(()=>en())}),En(()=>{u.value=window.innerWidth,window.addEventListener("resize",Ms),window.addEventListener("scroll",$s,{passive:!0})}),zn(()=>{window.removeEventListener("resize",Ms),window.removeEventListener("scroll",$s),h&&clearTimeout(h)}),(os,L)=>{const W=Za("router-link");return I(),N("div",c1,[S("div",{class:"reading-progress",style:qn({width:r.value+"%"}),"aria-hidden":"true"},null,4),S("div",r1,[S("div",i1,[S("div",{class:"sticky-sidebar",ref_key:"leftSidebarContent",ref:f},[m.value?(I(),N("div",u1,[bs(uu)])):ms("",!0)],512)]),S("div",h1,[S("div",d1,[S("div",m1,[S("div",j1,[P.value?(I(),N("div",g1,[S("h1",f1,G(P.value.title),1),S("div",b1,[g.value&&P.value.date?(I(),N("span",_1,[L[0]||(L[0]=S("i",{class:"fas fa-calendar-alt"},null,-1)),Cs(G(j.value)+" "+G(P.value.date),1)])):ms("",!0),Q.value>0?(I(),N("span",y1,[L[1]||(L[1]=S("i",{class:"fas fa-clock"},null,-1)),Cs(G(cs(Q.value)),1)])):ms("",!0),S("button",{type:"button",class:Ts(["article-copy-btn",{"article-copy-btn--copied":i.value}]),onClick:us,"aria-label":ss(n)("copyArticle"),"aria-live":"polite"},[S("i",{class:Ts(i.value?"fas fa-check":"fas fa-copy")},null,2),S("span",null,G(i.value?ss(n)("copied"):ss(n)("copyArticle")),1)],10,v1)]),P.value.tags?.length?(I(),N("div",w1,[(I(!0),N(gs,null,As(P.value.tags,(H,J)=>(I(),N("span",{key:J,class:"article-tag",onClick:hs=>X(H)}," # "+G(H),9,C1))),128))])):ms("",!0)])):ms("",!0),l.value?(I(),Aa(Uw,{key:1,rawMarkdown:l.value,articlePath:P.value?.path||"",articleTitle:P.value?.title||"",onMarkdownRendered:Mn},null,8,["rawMarkdown","articlePath","articleTitle"])):ms("",!0),l.value?(I(),N("nav",k1,[!V.value&&z.value?(I(),N("span",x1)):ms("",!0),V.value?(I(),Aa(W,{key:1,to:ls(V.value.path),class:"article-nav-item prev"},{default:jn(()=>[L[2]||(L[2]=S("div",{class:"nav-arrow"},[S("i",{class:"fas fa-arrow-left"})],-1)),S("div",S1,[S("div",P1,G(y.value),1),S("div",E1,G(V.value.title),1)])]),_:1},8,["to"])):ms("",!0),z.value?(I(),Aa(W,{key:2,to:ls(z.value.path),class:"article-nav-item next"},{default:jn(()=>[S("div",T1,[S("div",A1,G(x.value),1),S("div",R1,G(z.value.title),1)]),L[3]||(L[3]=S("div",{class:"nav-arrow"},[S("i",{class:"fas fa-arrow-right"})],-1))]),_:1},8,["to"])):ms("",!0)])):ms("",!0)])])])],512),S("div",L1,[S("div",{class:"sticky-sidebar",ref_key:"rightSidebarContent",ref:C},[m.value?(I(),N("div",M1,[bs(ju,{ref_key:"onThisPageRef",ref:b,containerSelector:".markdown-body",levels:[2,3],offset:88},null,512)])):ms("",!0)],512)])]),l.value?(I(),Aa(p1,{key:0})):ms("",!0)])}}}),O1=hn(D1,[["__scopeId","data-v-32d6ea44"]]),I1={class:"page-section notfound-view"},N1={class:"notfound-desc"},F1=Ns({name:"NotFoundView",head(){return{title:"404 - zorrooz’s blog"}},__name:"NotFound",setup(s){const{t:n}=Us(),a=ya(),e=ns(()=>`/${a.params.locale==="en"?"en":"zh"}/`);return(t,l)=>{const o=Za("router-link");return I(),N("div",I1,[l[1]||(l[1]=S("h1",{class:"notfound-code"},"404",-1)),S("p",N1,G(ss(n)("pageNotFound")),1),bs(o,{to:e.value,class:"notfound-back"},{default:jn(()=>[l[0]||(l[0]=S("i",{class:"fas fa-arrow-left"},null,-1)),Cs(G(ss(n)("backHome")),1)]),_:1},8,["to"])])}}}),$1=hn(F1,[["__scopeId","data-v-8d4bf42e"]]),B1=s=>[{path:`/${s}/`,name:`${s}-Home`,component:gv},{path:`/${s}/category`,name:`${s}-Category`,component:zv},{path:`/${s}/resource`,name:`${s}-Resource`,component:hw},{path:`/${s}/about`,name:`${s}-About`,component:qw},{path:`/${s}/article/:path*`,name:`${s}-Article`,component:O1,props:!0}],q1=[{path:"/",redirect:()=>`/${$a()}/`},{path:"/zh",redirect:"/zh/"},{path:"/en",redirect:"/en/"},{path:"/category",redirect:s=>`/${$a()}${s.path}`},{path:"/resource",redirect:s=>`/${$a()}${s.path}`},{path:"/about",redirect:s=>`/${$a()}${s.path}`},{path:"/article/:path*",redirect:s=>`/${$a()}${s.path}`},{path:"/:pathMatch(.*)*",redirect:s=>`/${$a()}${s.path}`},...pb.map(B1).flat(),{path:"/:locale(zh|en)/:pathMatch(.*)*",name:"NotFound",component:$1}];fs(()=>import("./bootstrap.esm-inoyDHN5.js"),[]);const z1={mounted(s){if(typeof window>"u"||window.matchMedia("(prefers-reduced-motion: reduce)").matches||!("IntersectionObserver"in window))return;s.classList.add("reveal");const n=new IntersectionObserver(([a])=>{a.isIntersecting&&(s.classList.add("reveal-visible"),n.disconnect())},{threshold:.12});n.observe(s),s.__revealIo=n},unmounted(s){s.__revealIo?.disconnect(),s.__revealIo=null}};Mj(Cy,{routes:q1},({app:s,router:n,initialState:a,isClient:e,routePath:t})=>{const l=Wd();s.use(l),a.pinia&&(l.state.value=a.pinia),s.use(n),s.use(Le),s.mixin(Cm),s.directive("reveal",z1);const o=Do(),p=Mo();o.initTheme(),e&&p.initLocale(),e&&"serviceWorker"in navigator&&window.addEventListener("load",()=>{navigator.serviceWorker.register("/sw.js").catch(()=>{})})});export{gs as F,V1 as T,fs as _,qe as a,S as b,Aa as c,Ns as d,ss as e,an as f,N as g,Cs as h,ms as i,As as j,Sp as k,vn as l,Jn as m,bn as n,I as o,hn as p,rs as r,G as t,Us as u,U1 as v,Vs as w};
