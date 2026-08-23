<template>
  <nav class="on-this-page">
    <div class="otp-content" ref="otpContent">
      <span class="otp-marker" ref="marker" aria-hidden="true"></span>
      <div class="otp-header">
        <span class="otp-title">{{ t('tableOfContents') }}</span>
      </div>

      <ul class="otp-list">
        <li
          v-for="(item, idx) in toc"
          :key="idx"
          :class="['otp-item', { active: activeId === item.id }]"
        >
          <a
            class="otp-link"
            role="button"
            tabindex="0"
            @click.prevent="scrollToId(item.id)"
            @keydown.enter.prevent="scrollToId(item.id)"
          >
            <span class="otp-text">{{ item.text }}</span>
          </a>
          <ul v-if="item.children && item.children.length" class="otp-sublist">
            <li
              v-for="(sub, sIdx) in item.children"
              :key="sIdx"
              :class="['otp-subitem', { active: activeId === sub.id }]"
            >
              <a
                class="otp-sublink"
                role="button"
                tabindex="0"
                @click.prevent="scrollToId(sub.id)"
                @keydown.enter.prevent="scrollToId(sub.id)"
              >
                <span class="otp-subtext">{{ sub.text }}</span>
              </a>
            </li>
          </ul>
        </li>
      </ul>
    </div>
  </nav>
</template>

<script setup lang="ts">
import { nextTick, onBeforeUnmount, onMounted, ref, useTemplateRef, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { scrollToHeading } from '@/utils/scroll'
import { buildTocTree, ensureHeadingIds, type TocNode } from '@/utils/toc'

const props = withDefaults(
  defineProps<{
    containerSelector?: string
    levels?: number[]
    offset?: number
  }>(),
  {
    containerSelector: '.markdown-body',
    levels: () => [2, 3, 4, 5, 6],
    offset: 8,
  },
)

const emit = defineEmits(['navigate'])

const { t } = useI18n()

const toc = ref<TocNode[]>([])
const activeId = ref('')
const otpContent = useTemplateRef<HTMLElement>('otpContent')
const marker = useTemplateRef<HTMLElement>('marker')

/** 标题行高 32px + marker 相对链接首行的居中偏移，与 VitePress outline-marker 一致 */
const MARKER_TITLE_OFFSET = 36

/** observer/poller 为非响应式句柄，普通变量持有即可 */
let containerObserver: MutationObserver | null = null
let containerTimer: number | null = null
let containerPoller: number | null = null
let spyTicking = false

/**
 * 点击导航后的抑制窗口：平滑滚动途中高亮会随经过的标题短暂变化，
 * 窗口期内以点击目标为准，结束后 scrollspy 自动恢复接管。
 */
const SPY_SUPPRESS_MS = 900
let spySuppressedUntil = 0

function cleanupObservers() {
  containerObserver?.disconnect()
  containerObserver = null
  if (containerTimer) {
    clearTimeout(containerTimer)
    containerTimer = null
  }
  if (containerPoller) {
    clearInterval(containerPoller)
    containerPoller = null
  }
  window.removeEventListener('scroll', onSpyScroll)
  window.removeEventListener('resize', onSpyScroll)
}

/**
 * 滚动位置驱动的 scrollspy：取「判定线（视口顶部 + offset）以上最近的标题」。
 * 不用 IntersectionObserver——目标标题被 scrollToHeading 定位在观察区上边缘时
 * 会被判为不相交，高亮被其下第一个子标题抢走。
 */
function updateActiveByScroll() {
  if (Date.now() < spySuppressedUntil) return
  const root = document.querySelector(props.containerSelector)
  if (!root) return
  const selector = props.levels.map((l) => `h${l}`).join(',')
  const headings = root.querySelectorAll(selector)
  if (headings.length === 0) return

  const line = window.scrollY + props.offset + 1
  let current = ''
  for (const h of headings) {
    if (h.getBoundingClientRect().top + window.scrollY <= line) current = h.id
    else break
  }
  activeId.value = current
}

function onSpyScroll() {
  if (spyTicking) return
  spyTicking = true
  requestAnimationFrame(() => {
    spyTicking = false
    updateActiveByScroll()
  })
}

function setupScrollSpy() {
  updateActiveByScroll()
}

function resetToc() {
  toc.value = []
  activeId.value = ''
  cleanupObservers()
  nextTick(() => setupContainerObserver())
}

function refreshToc() {
  buildToc()
  setupScrollSpy()
  nextTick(updateMarker)
}

/** marker 滑动到当前激活项（offsetTop 相对 .otp-list，+标题行高得到容器内位置） */
function updateMarker() {
  if (!marker.value || !otpContent.value) return
  const root = otpContent.value.querySelector('.otp-list')
  if (!root) return
  const activeLink = root.querySelector<HTMLElement>(
    '.otp-item.active > .otp-link, .otp-subitem.active > .otp-sublink',
  )
  if (!activeLink) {
    marker.value.style.opacity = '0'
    return
  }
  marker.value.style.top = `${activeLink.offsetTop + MARKER_TITLE_OFFSET}px`
  marker.value.style.opacity = '1'
}

watch(activeId, () => nextTick(updateMarker))

/** 容器晚于组件挂载时轮询等待，就绪后以 MutationObserver 跟随内容变化 */
function setupContainerObserver() {
  setTimeout(() => refreshToc(), 100)
  const checkContainer = () => {
    const root = document.querySelector(props.containerSelector)
    if (!root) return
    if (containerPoller) {
      clearInterval(containerPoller)
      containerPoller = null
    }
    if (!containerObserver) {
      containerObserver = new MutationObserver(() => {
        clearTimeout(containerTimer ?? undefined)
        containerTimer = window.setTimeout(() => refreshToc(), 100)
      })
      containerObserver.observe(root, {
        childList: true,
        subtree: true,
        attributes: true,
        attributeFilter: ['id'],
      })
    }
    refreshToc()
  }
  checkContainer()
  if (!containerPoller && !document.querySelector(props.containerSelector))
    containerPoller = window.setInterval(checkContainer, 200)
}

function buildToc() {
  const root = document.querySelector(props.containerSelector)
  if (!root) {
    toc.value = []
    return
  }

  const selector = props.levels.map((l) => `h${l}`).join(',')
  const headings = Array.from(root.querySelectorAll(selector))
  if (headings.length === 0) {
    toc.value = []
    return
  }

  ensureHeadingIds(headings)
  toc.value = buildTocTree(headings, props.levels)
}

function scrollToId(id: string) {
  emit('navigate', id)

  const el = document.getElementById(id)
  if (!el) return
  // 点击目标立即高亮；抑制窗口过后 scrollspy 自动恢复接管
  activeId.value = id
  spySuppressedUntil = Date.now() + SPY_SUPPRESS_MS
  scrollToHeading(el, props.offset)
}

onMounted(() => {
  buildToc()
  setupScrollSpy()
  window.addEventListener('scroll', onSpyScroll, { passive: true })
  window.addEventListener('resize', onSpyScroll)
  nextTick(() => {
    updateMarker()
    setupContainerObserver()
  })
})

onBeforeUnmount(() => {
  cleanupObservers()
})

defineExpose({ refreshToc, resetToc })
</script>

<style scoped>
.on-this-page {
  padding: 0.25rem 0;
  font-size: var(--text-base);
}

/* VitePress outline 容器：左侧发丝导轨 + 滑动 marker */
.otp-content {
  position: relative;
  border-left: 1px solid var(--line);
  padding-left: 16px;
  font-size: 13px;
  font-weight: 500;
}

.otp-marker {
  position: absolute;
  top: 32px;
  left: -1px;
  z-index: 0;
  width: 2px;
  height: 18px;
  border-radius: 2px;
  background: var(--primary);
  opacity: 0;
  pointer-events: none;
  transition:
    top 0.25s cubic-bezier(0, 1, 0.5, 1),
    opacity 0.25s ease;
}

.otp-header {
  line-height: 32px;
}

.otp-title {
  font-size: 13px;
  font-weight: 600;
  color: var(--fg);
}

.otp-list {
  list-style: none;
  padding: 0;
  margin: 0;
  position: relative;
  z-index: 1;
}

.otp-item {
  margin: 0;
}

/* 换行显示（不截断省略），多行时 marker 锚定首行 */
.otp-link {
  display: block;
  line-height: 1.55;
  padding: 3px 0;
  font-size: 13px;
  font-weight: 400;
  color: var(--fg-2);
  text-decoration: none;
  transition: color var(--dur-base) ease;
  cursor: pointer;
}

.otp-link:hover,
.otp-item.active > .otp-link {
  color: var(--fg);
  transition: color var(--dur-fast) ease;
}

.otp-sublist {
  list-style: none;
  padding: 0;
  margin: 0;
}

.otp-sublink {
  display: block;
  padding: 3px 0 3px 13px;
  line-height: 1.55;
  font-size: 13px;
  font-weight: 400;
  color: var(--fg-2);
  text-decoration: none;
  transition: color var(--dur-base) ease;
  cursor: pointer;
}

.otp-sublink:hover,
.otp-subitem.active > .otp-sublink {
  color: var(--fg);
  transition: color var(--dur-fast) ease;
}
</style>
