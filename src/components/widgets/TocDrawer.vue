<!-- TocDrawer.vue -->
<template>
  <button v-show="visible" class="toc-drawer-btn d-lg-none d-flex align-items-center justify-content-center"
    @click="openDrawer" :aria-label="t('openToc')" @touchstart.prevent.stop="onTouchStart"
    @touchmove.prevent.stop="onTouchMove" @touchend.prevent.stop="onTouchEnd" :style="{ top: buttonTop + 'px' }">
    <i class="fas fa-bookmark"></i>
  </button>

  <div v-if="show" class="mobile-offcanvas d-lg-none" @click.self="close">
    <div class="offcanvas-panel offcanvas-right border-start rounded-0 shadow-sm">
      <div class="offcanvas-section">
        <div class="offcanvas-card">
          <OnThisPage containerSelector=".markdown-body" :levels="[2, 3]" :offset="8" @navigate="onNavigate" />
        </div>
      </div>
    </div>
    <div class="offcanvas-backdrop" @click="close"></div>
  </div>
</template>

<script>
/*
  TocDrawer
  - 移动端目录抽屉，右侧面板，支持拖动调整按钮位置并与 BackToTop 同步
*/
import { useI18n } from 'vue-i18n'
import OnThisPage from '@/components/layout/OnThisPage.vue';

export default {
  name: 'TocDrawer',
  setup() { const { t } = useI18n(); return { t } },
  components: { OnThisPage },
  data() {
    return {
      sourceId: 'toc', rafPending: false, rafLastBaseTop: null,
      visible: false, isDragging: false, startY: 0, initialTop: 0, buttonTop: window.innerHeight - 160, touchMoved: false,
      show: false,
      _lockedScrollY: null
    };
  },
  mounted() {
    window.addEventListener('scroll', this.onScroll, { passive: true });
    window.addEventListener('resize', this.onResize, { passive: true });
    window.addEventListener('floating-buttons-base-top', this.syncBaseTop);
    this.onResize(); this.onScroll();
    const GAP = 48; this.rafDispatchBaseTop(this.buttonTop + GAP);
  },
  beforeUnmount() {
    window.removeEventListener('scroll', this.onScroll);
    window.removeEventListener('resize', this.onResize);
    window.removeEventListener('floating-buttons-base-top', this.syncBaseTop);
    this.unlockScroll();
  },
  methods: {
    getBounds() { const BUTTON_HEIGHT = 40; const GAP = 48; const MARGIN = 20; return { gap: GAP, minTop: MARGIN, maxTop: Math.max(this.$el ? 0 : 0, window.innerHeight - BUTTON_HEIGHT - MARGIN - GAP) } },
    clampTop(top) { const { minTop, maxTop } = this.getBounds(); return Math.max(minTop, Math.min(maxTop, top)); },
    rafDispatchBaseTop(baseTop) { this.rafLastBaseTop = baseTop; if (this.rafPending) return; this.rafPending = true; requestAnimationFrame(() => { window.dispatchEvent(new CustomEvent('floating-buttons-base-top', { detail: { baseTop: this.rafLastBaseTop, source: this.sourceId } })); this.rafPending = false; }); },
    syncBaseTop(e) { const base = e?.detail?.baseTop; const source = e?.detail?.source; if (source === this.sourceId) return; const { gap } = this.getBounds(); if (typeof base === 'number') { const desiredTop = base - gap; this.buttonTop = this.clampTop(desiredTop); } },
    onResize() { const isMobile = window.innerWidth < 992; this.visible = isMobile; this.buttonTop = this.clampTop(this.buttonTop); },
    onScroll() { if (!this.isDragging) { const isMobile = window.innerWidth < 992; this.visible = isMobile; } },
    openDrawer() { if (this.isDragging) return; this.lockScroll(); this.show = true; },
    close() { this.show = false; this.unlockScroll(); },
    onTouchStart(e) { e.preventDefault(); this.isDragging = true; this.touchMoved = false; this.startY = e.touches[0].clientY; this.initialTop = this.buttonTop; },
    onTouchMove(e) { this.touchMoved = true; const currentY = e.touches[0].clientY; const diffY = currentY - this.startY; let newTop = this.clampTop(this.initialTop + diffY); this.buttonTop = newTop; const { gap } = this.getBounds(); this.rafDispatchBaseTop(newTop + gap); e.preventDefault(); },
    onTouchEnd(e) { e.preventDefault(); this.isDragging = false; if (!this.touchMoved) this.openDrawer(); },
    onNavigate(id) { this.close(); },
    lockScroll() {
      try {
        const scrollY = window.scrollY || window.pageYOffset || document.documentElement.scrollTop || 0;
        this._lockedScrollY = scrollY;
        const body = document.body;
        if (body) {
          body.style.position = 'fixed';
          body.style.top = `-${scrollY}px`;
          body.style.left = '0';
          body.style.right = '0';
          body.style.overflow = 'hidden';
        }
        if (document.documentElement) {
          document.documentElement.style.overscrollBehavior = 'contain';
        }
      } catch (e) {
        const docEl = document.documentElement;
        const body = document.body;
        if (docEl) { docEl.style.overflow = 'hidden'; docEl.style.overscrollBehavior = 'contain'; }
        if (body) { body.style.overflow = 'hidden'; body.style.overscrollBehavior = 'contain'; }
      }
    },
    unlockScroll() {
      try {
        const body = document.body;
        if (body) {
          body.style.position = '';
          body.style.top = '';
          body.style.left = '';
          body.style.right = '';
          body.style.overflow = '';
        }
        if (document.documentElement) {
          document.documentElement.style.overscrollBehavior = '';
        }
        if (typeof this._lockedScrollY === 'number') {
          window.scrollTo(0, this._lockedScrollY);
          this._lockedScrollY = null;
        }
      } catch (e) {
        const docEl = document.documentElement;
        const body = document.body;
        if (docEl) { docEl.style.overflow = ''; docEl.style.overscrollBehavior = ''; }
        if (body) { body.style.overflow = ''; body.style.overscrollBehavior = ''; }
      }
    }
  }
};
</script>

<style scoped>
.toc-drawer-btn {
  position: fixed;
  right: 30px;
  width: 40px;
  height: 40px;
  background-color: var(--app-btn-bg);
  color: var(--app-text);
  border: 1px solid var(--app-border);
  border-radius: 8px;
  font-size: 18px;
  font-weight: bold;
  cursor: pointer;
  box-shadow: var(--app-btn-shadow);
  z-index: 1000;
  outline: none;
  -webkit-tap-highlight-color: transparent;
  touch-action: none;
}

.toc-drawer-btn:hover {
  background-color: var(--app-btn-hover-bg);
  box-shadow: var(--app-btn-hover-shadow);
}

.toc-drawer-btn:active {
  transform: scale(0.95);
}

.mobile-offcanvas {
  position: fixed;
  inset: 0;
  z-index: 1050;
}

.offcanvas-backdrop {
  position: absolute;
  inset: 0;
  background: var(--app-backdrop-bg);
  z-index: 1;
}

.offcanvas-panel {
  position: absolute;
  top: 0;
  bottom: 0;
  width: min(85vw, 320px);
  background: var(--app-offcanvas-bg);
  box-shadow: var(--app-offcanvas-shadow);
  overflow-y: auto;
  -webkit-overflow-scrolling: touch;
  padding: 0.75rem 0.75rem 1rem;
  z-index: 2;
}

.offcanvas-right {
  right: 0;
}

.offcanvas-section+.offcanvas-section {
  margin-top: 0.75rem;
}

.section-title {
  font-size: 0.95rem;
  color: var(--app-text-muted);
  margin-bottom: 0.25rem;
  font-weight: 600;
}

.offcanvas-card {
  background-color: var(--app-card-bg);
  border: none;
  border-radius: 0.5rem;
  padding: 0.5rem;
  margin: 1rem 0;
}
</style>