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
          <OnThisPage containerSelector=".markdown-body" :levels="[2, 3]" :offset="88" @navigate="onNavigate" />
        </div>
      </div>
    </div>
    <div class="offcanvas-backdrop" @click="close"></div>
  </div>
</template>

<script setup lang="ts">
import { onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import OnThisPage from '@/components/layout/OnThisPage.vue'

const { t } = useI18n()

const sourceId = 'toc'
const rafPending = ref(false)
const rafLastBaseTop = ref<number | null>(null)
const visible = ref(false)
const isDragging = ref(false)
const startY = ref(0)
const initialTop = ref(0)
const buttonTop = ref((typeof window !== 'undefined' ? window.innerHeight : 1024) - 160)
const touchMoved = ref(false)
const show = ref(false)
const lockedScrollY = ref<number | null>(null)

function getBounds() {
  const BUTTON_HEIGHT = 40
  const GAP = 48
  const MARGIN = 20
  return { gap: GAP, minTop: MARGIN, maxTop: Math.max(0, window.innerHeight - BUTTON_HEIGHT - MARGIN - GAP) }
}

function clampTop(top: number) {
  const { minTop, maxTop } = getBounds()
  return Math.max(minTop, Math.min(maxTop, top))
}

function rafDispatchBaseTop(baseTop: number) {
  rafLastBaseTop.value = baseTop
  if (rafPending.value) return
  rafPending.value = true
  requestAnimationFrame(() => {
    window.dispatchEvent(new CustomEvent('floating-buttons-base-top', { detail: { baseTop: rafLastBaseTop.value, source: sourceId } }))
    rafPending.value = false
  })
}

function syncBaseTop(e: Event) {
  const detail = (e as CustomEvent).detail
  const base = detail?.baseTop
  const source = detail?.source
  if (source === sourceId) return
  const { gap } = getBounds()
  if (typeof base === 'number') {
    const desiredTop = base - gap
    buttonTop.value = clampTop(desiredTop)
  }
}

function onResize() {
  const isMobile = window.innerWidth < 992
  visible.value = isMobile
  buttonTop.value = clampTop(buttonTop.value)
}

function openDrawer() {
  if (isDragging.value) return
  lockScroll()
  show.value = true
}

function close() {
  show.value = false
  unlockScroll()
}

function onTouchStart(e: TouchEvent) {
  e.preventDefault(); isDragging.value = true; touchMoved.value = false; startY.value = e.touches[0].clientY; initialTop.value = buttonTop.value
}

function onTouchMove(e: TouchEvent) {
  touchMoved.value = true
  const currentY = e.touches[0].clientY
  const diffY = currentY - startY.value
  const newTop = clampTop(initialTop.value + diffY)
  buttonTop.value = newTop
  const { gap } = getBounds()
  rafDispatchBaseTop(newTop + gap)
  e.preventDefault()
}

function onTouchEnd(e: TouchEvent) {
  e.preventDefault(); isDragging.value = false; if (!touchMoved.value) openDrawer()
}

function onNavigate() { close() }

function lockScroll() {
  try {
    const scrollY = window.scrollY || window.pageYOffset || document.documentElement.scrollTop || 0
    lockedScrollY.value = scrollY
    const body = document.body
    if (body) {
      body.style.position = 'fixed'
      body.style.top = `-${scrollY}px`
      body.style.left = '0'
      body.style.right = '0'
      body.style.overflow = 'hidden'
    }
    if (document.documentElement) {
      document.documentElement.style.overscrollBehavior = 'contain'
    }
  } catch {
    const docEl = document.documentElement
    const body = document.body
    if (docEl) { docEl.style.overflow = 'hidden'; docEl.style.overscrollBehavior = 'contain' }
    if (body) { body.style.overflow = 'hidden'; body.style.overscrollBehavior = 'contain' }
  }
}

function unlockScroll() {
  try {
    const body = document.body
    if (body) {
      body.style.position = ''
      body.style.top = ''
      body.style.left = ''
      body.style.right = ''
      body.style.overflow = ''
    }
    if (document.documentElement) {
      document.documentElement.style.overscrollBehavior = ''
    }
    if (typeof lockedScrollY.value === 'number') {
      window.scrollTo(0, lockedScrollY.value)
      lockedScrollY.value = null
    }
  } catch {
    const docEl = document.documentElement
    const body = document.body
    if (docEl) { docEl.style.overflow = ''; docEl.style.overscrollBehavior = '' }
    if (body) { body.style.overflow = ''; body.style.overscrollBehavior = '' }
  }
}

onMounted(() => {
  window.addEventListener('resize', onResize, { passive: true })
  window.addEventListener('floating-buttons-base-top', syncBaseTop)
  onResize()
  const GAP = 48
  rafDispatchBaseTop(buttonTop.value + GAP)
})

onBeforeUnmount(() => {
  window.removeEventListener('resize', onResize)
  window.removeEventListener('floating-buttons-base-top', syncBaseTop)
  unlockScroll()
})
</script>

<style scoped>
.toc-drawer-btn {
  position: fixed;
  right: 28px;
  width: 40px;
  height: 40px;
  background-color: var(--surface);
  color: var(--fg-2);
  border: 1px solid var(--line);
  border-radius: 50%;
  font-size: 14px;
  cursor: pointer;
  box-shadow: var(--shadow-soft);
  z-index: 1000;
  outline: none;
  -webkit-tap-highlight-color: transparent;
  touch-action: none;
  transition: background-color 0.14s ease, box-shadow 0.14s ease, transform 0.14s ease,
    color 0.14s ease, border-color 0.14s ease;
}

.toc-drawer-btn:hover {
  background-color: var(--primary);
  color: var(--on-primary);
  border-color: transparent;
  box-shadow: var(--shadow-lift);
}

.toc-drawer-btn:active {
  transform: scale(0.93);
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
  width: min(82vw, 300px);
  background: var(--surface);
  border-left: 1px solid var(--line);
  box-shadow: var(--shadow-lift);
  overflow-y: auto;
  -webkit-overflow-scrolling: touch;
  padding: 1rem 0.75rem 1.5rem;
  z-index: 2;
}

.offcanvas-right {
  right: 0;
}

.offcanvas-section + .offcanvas-section {
  margin-top: 0.75rem;
}

.offcanvas-card {
  background-color: var(--surface);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  padding: 0.5rem;
  margin: 0.5rem 0;
}
</style>
