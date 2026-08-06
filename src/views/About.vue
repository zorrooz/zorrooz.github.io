<template>
  <div class="page-section about-view">
    <header class="about-head">
      <div class="about-head__identity">
        <div class="about-head__avatar">
          <img
            v-if="avatarSrc"
            :src="avatarSrc"
            alt="avatar"
            class="about-head__avatar-img"
            draggable="false"
          />
          <span v-else class="about-head__initial">{{ siteAuthorInitial }}</span>
        </div>
        <div class="about-head__names">
          <h1 class="about-head__name">{{ siteAuthor }}</h1>
          <p class="about-head__role">{{ t('developer') }}</p>
        </div>
      </div>

      <p v-if="introduction" class="about-intro">{{ introduction }}</p>
    </header>

    <main class="about-body">
      <section v-if="experience.length" class="about-section">
        <div class="about-section__title">{{ t('experience') }}</div>
        <div class="timeline">
          <div v-for="(it, idx) in experience" :key="idx" class="tl-item">
            <div class="tl-year num">{{ it.year }}</div>
            <div class="tl-content">
              <div class="tl-title">{{ it.title }}</div>
              <div v-if="it.desc" class="tl-desc">{{ it.desc }}</div>
            </div>
          </div>
        </div>
      </section>

      <div class="about-grid" v-if="sections.length">
        <div v-for="(sec, idx) in sections" :key="idx" class="about-cell">
          <div class="about-cell__title">{{ sec.title }}</div>
          <div v-for="(it, j) in sec.items" :key="j" class="about-cell__item">
            <span v-if="it.name" class="about-cell__item-name">{{ it.name }}</span>
            <span v-if="it.desc" class="about-cell__item-desc">{{ it.desc }}</span>
          </div>
        </div>
      </div>
    </main>

    <footer class="about-foot">
      <div class="about-foot__contacts">
        <a
          v-for="contact in data.contacts"
          :key="contact.label"
          :href="contact.link"
          target="_blank"
          rel="noopener noreferrer"
          class="about-foot__contact"
        >
          <i :class="contact.icon"></i>
          <span>{{ contact.value }}</span>
        </a>
      </div>
      <p class="about-foot__thanks">{{ footerText }}</p>
    </footer>
  </div>
</template>

<script setup lang="ts">
defineOptions({ name: 'AboutView' })
import { computed } from 'vue'
import { useI18n } from 'vue-i18n'
import { useLocalizedContent } from '@/composables/useLocalizedContent'
import { EMPTY_ABOUT, loadAbout } from '@/utils/contentLoader'
import { SITE } from '@/config'

const { t } = useI18n()

const { data } = useLocalizedContent(() => loadAbout(), EMPTY_ABOUT)

const avatarModules = import.meta.glob('/src/assets/avatar.*', {
  query: '?url',
  import: 'default',
  eager: true,
})
const avatarSrc = computed(() => (Object.values(avatarModules)[0] as string | undefined) || '')

const footerText = computed(() => t('thanks'))
const siteAuthor = SITE.author
const siteAuthorInitial = siteAuthor.trim().charAt(0).toUpperCase()

const introduction = computed(() => data.value.introduction)
const experience = computed(() => data.value.experience)
const sections = computed(() => data.value.section)
</script>

<style scoped>
/* ---------- Head: centered identity + intro ---------- */
.about-head {
  display: flex;
  flex-direction: column;
  align-items: center;
  text-align: center;
  gap: var(--sp-8);
  padding: var(--sp-16) 0 var(--sp-16);
  border-bottom: 1px solid var(--line);
  margin-bottom: var(--sp-16);
}

.about-head__identity {
  display: flex;
  align-items: center;
  justify-content: center;
  gap: var(--sp-6);
}

.about-head__avatar {
  width: 88px;
  height: 88px;
  border-radius: var(--radius-lg);
  background: transparent;
  border: 1.5px solid var(--ink);
  display: flex;
  align-items: center;
  justify-content: center;
  color: var(--ink);
  flex-shrink: 0;
  overflow: hidden;
}

.about-head__avatar-img {
  width: 100%;
  height: 100%;
  object-fit: cover;
  display: block;
}

.about-head__initial {
  font-family: var(--font-serif);
  font-size: 40px;
  font-weight: 700;
  line-height: 1;
}

.about-head__name {
  font-size: clamp(32px, 5vw, 48px);
  font-weight: 700;
  letter-spacing: -0.03em;
  color: var(--fg);
  margin: 0;
  line-height: 1.1;
}

.about-head__role {
  font-size: var(--text-sm);
  font-weight: 600;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--fg-3);
  margin: var(--sp-2) 0 0;
}

.about-intro {
  font-size: var(--text-body);
  color: var(--fg-2);
  line-height: 1.75;
  margin: 0 auto;
  max-width: 65ch;
}

/* ---------- Foot: contacts ---------- */
.about-foot {
  display: flex;
  flex-direction: column;
  align-items: center;
  text-align: center;
  gap: var(--sp-4);
  margin-top: var(--sp-16);
  padding-top: var(--sp-8);
  border-top: 1px solid var(--line);
}

.about-foot__contacts {
  display: flex;
  flex-wrap: wrap;
  justify-content: center;
  gap: var(--sp-2) var(--sp-8);
}

.about-foot__contact {
  display: inline-flex;
  align-items: center;
  gap: 8px;
  color: var(--fg-2);
  font-size: var(--text-sm);
  font-weight: 500;
  text-decoration: none;
  transition: color 0.14s ease;
}

.about-foot__contact i {
  color: var(--fg-3);
  font-size: 14px;
  transition: color 0.14s ease;
}

.about-foot__contact:hover {
  color: var(--primary);
}

.about-foot__contact:hover i {
  color: var(--primary);
}

.about-foot__thanks {
  margin: 0;
  font-size: var(--text-xs);
  color: var(--fg-3);
}

/* ---------- Body ---------- */
.about-body {
  display: flex;
  flex-direction: column;
  gap: var(--sp-16);
}

/* ---------- Experience (centered) ---------- */
.about-section {
  width: 100%;
  max-width: 720px;
  margin-inline: auto;
}

.about-section__title {
  display: flex;
  align-items: center;
  gap: var(--sp-3);
  font-size: var(--text-base);
  font-weight: 600;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--fg-3);
  margin-bottom: var(--sp-6);
}

.about-section__title::before {
  content: '';
  width: 3px;
  height: 16px;
  border-radius: 99px;
  background: var(--primary);
  flex-shrink: 0;
}

/* ---------- Timeline ---------- */
.timeline {
  padding: var(--sp-2) 0;
}

.tl-item {
  position: relative;
  padding: 0 0 var(--sp-8) var(--sp-6);
  border-left: 1px solid var(--line);
  transition: background-color 0.14s ease;
}

.tl-item:last-child {
  padding-bottom: var(--sp-3);
}

.tl-item::before {
  content: '';
  position: absolute;
  left: -4px;
  top: 5px;
  width: 7px;
  height: 7px;
  border-radius: 50%;
  background: var(--primary);
  transition:
    transform 0.14s ease,
    box-shadow 0.14s ease;
}

.tl-item:hover {
  background: var(--tint);
  border-radius: 0 var(--radius) var(--radius) 0;
}

.tl-item:hover::before {
  transform: scale(1.35);
  box-shadow: 0 0 0 4px color-mix(in srgb, var(--primary) 15%, transparent);
}

.tl-year {
  font-family: 'Agave', 'SourceHanSansSC', var(--font-sans);
  font-size: var(--text-sm);
  font-weight: 600;
  letter-spacing: 0.04em;
  color: var(--primary);
  margin-bottom: var(--sp-2);
  line-height: 1.5;
  font-variant-numeric: tabular-nums;
}

.tl-title {
  font-size: 18px;
  font-weight: 600;
  color: var(--fg);
  line-height: 1.5;
  transition: color 0.14s ease;
}

.tl-item:hover .tl-title {
  color: var(--primary);
}

.tl-desc {
  font-size: var(--text-base);
  color: var(--fg-2);
  margin-top: 3px;
  line-height: 1.65;
}

/* ---------- Section cards ---------- */
.about-grid {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: var(--sp-5);
  align-content: start;
}

.about-cell {
  background: var(--surface);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  padding: var(--sp-8);
  transition:
    border-color 0.18s ease,
    box-shadow 0.18s ease;
}

.about-cell:hover {
  border-color: color-mix(in srgb, var(--primary) 30%, transparent);
  box-shadow: var(--shadow-soft);
}

.about-cell__title {
  font-size: var(--text-sm);
  font-weight: 600;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--primary);
  margin-bottom: var(--sp-4);
  display: flex;
  align-items: center;
  gap: 8px;
}

.about-cell__title::before {
  content: '';
  width: 2px;
  height: 14px;
  border-radius: 99px;
  background: var(--primary);
  flex-shrink: 0;
}

.about-cell__item {
  display: flex;
  flex-direction: column;
  margin-bottom: var(--sp-3);
  line-height: 1.6;
}

.about-cell__item:last-child {
  margin-bottom: 0;
}

.about-cell__item-name {
  font-size: var(--text-md);
  font-weight: 600;
  color: var(--fg);
}

.about-cell__item-desc {
  font-size: var(--text-base);
  color: var(--fg-3);
  margin-top: 2px;
}

/* ---------- Responsive ---------- */
@media (max-width: 767px) {
  .about-grid {
    grid-template-columns: 1fr;
  }

  .about-body {
    gap: var(--sp-12);
  }

  .about-head {
    padding: var(--sp-12) 0 var(--sp-12);
    gap: var(--sp-6);
  }

  .about-head__avatar {
    width: 72px;
    height: 72px;
  }

  .about-head__initial {
    font-size: 34px;
  }

  .about-head__name {
    font-size: clamp(28px, 8vw, 36px);
  }

  .about-foot__contacts {
    gap: var(--sp-3) var(--sp-6);
  }
}
</style>
